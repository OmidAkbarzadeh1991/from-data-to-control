clc
clear
close all
tic

ops.solver = 'mosek';

rng(1)

kappa = 0.001;
T     = 15;
tau   = 0.01;
Nvals = [10 50 100 500 1000];
n     = 3;

X_bound = [-20 20; -20 20; -20 20];

sigma = 10;
r     = 28;
b     = 8/3;

A = [-sigma, sigma,  0,  0, 0;
      r,    -1,      0, -1, 0;
      0,     0,     -b,  0, 1];
B = eye(3);
D = kappa * eye(3);

initial = [2; 2; 2];
states = initial.'; statesdot = []; N0T = []; u = []; w = [];
for i = 1:T
    x0  = initial;
    M   = [x0(1), x0(2), x0(3), x0(1)*x0(3), x0(1)*x0(2)];
    M_c = [x0(1), x0(2), x0(3), x0(1)*x0(3), x0(1)*x0(2), ...
           x0(2)*x0(3), x0(2)^2, x0(1)^2, x0(3)^2];
    u(i, :) = 100 * (-1 + 2*rand(1, 3));
    w(i, :) = 10 * (-1 + 2*rand(1, 3));
    odeSystem = @(t, x) A*M.' + B*u(i, :).' + D*w(i, :).';
    [~, X]    = ode45(odeSystem, [(i-1)*tau, i*tau], x0);
    initial   = X(end, :).';
    states    = [states; initial.'];
    statesdot = [statesdot; (A*M.' + B*u(i, :).' + D*w(i, :).').'];
    N0T       = [N0T; M_c];
end
U0T = u.'; W0T = w.'; X0T = states(1:end-1, :).';
X1T = statesdot.' + (-0.1 + 0.2*rand(3, T));
N0T = N0T.';
bar_varkappa = 3 * 0.01;

pvar x_1 x_2 x_3
x = [x_1; x_2; x_3];

Upsilon_x = [1 0 0; 0 1 0; 0 0 1; x_3 0 0; x_2 0 0; 0 x_3 0; ...
             0 x_2 0; x_1 0 0; 0 0 x_3];

prog = sosprogram(x);
Mon_Lagrangian = monomials(x, 0:2);
b_x = [(x_1 - X_bound(1, 1)) * (X_bound(1, 2) - x_1);
       (x_2 - X_bound(2, 1)) * (X_bound(2, 2) - x_2);
       (x_3 - X_bound(3, 1)) * (X_bound(3, 2) - x_3)];
[prog, lx1] = sospolyvar(prog, Mon_Lagrangian, 'wscoeff');
[prog, lx2] = sospolyvar(prog, Mon_Lagrangian, 'wscoeff');
[prog, lx3] = sospolyvar(prog, Mon_Lagrangian, 'wscoeff');
lambda_x = [lx1; lx2; lx3];

PI = 0.8481;
[prog, mu]   = sospolyvar(prog, monomials(x, 0));
[prog, invP] = sospolymatrixvar(prog, monomials(x, 0), [n n], 'symmetric');
[prog, H_x]  = sospolymatrixvar(prog, monomials(x, 0:1), [T n]);
epsilon = 0.99;

barx11 = (1/PI) * (D.' * D);
barx21 = zeros(n);
[prog, barx22] = sospolymatrixvar(prog, monomials(x, 0), [n n], 'symmetric');

G = - epsilon * invP ...
    - ((X1T - D*W0T) * H_x) - (H_x.' * (X1T - D*W0T).') ...
    - mu * (T * bar_varkappa * eye(n)) ...
    - PI * eye(n);

Condition_3 = [ G + barx22,  H_x.';
                H_x,         mu * eye(T) ] ...
              - lambda_x.' * b_x * eye(T + n);

for ii = 1:length(lambda_x), prog = sosineq(prog, lambda_x(ii)); end
prog = soseq(prog, N0T * H_x - Upsilon_x * invP);
prog = sosmatrixineq(prog, invP - 1e-6 * eye(n), 'quadraticMineq');
prog = sosmatrixineq(prog, Condition_3, 'quadraticMineq');
prog = sosmatrixineq(prog, -barx22 - eye(n), 'quadraticMineq');
prog = sosineq(prog, mu - 1e-6);
prog = sossolve(prog, ops);

info = prog.solinfo.info;
assert(info.pinf == 0 && info.dinf == 0 && info.numerr == 0, ...
       'Dissipativity synthesis infeasible; relax barx22 <= -I to e.g. -0.1*I.');

mu     = double(sosgetsol(prog, mu));
invP   = double(sosgetsol(prog, invP));
P      = invP^(-1);
H_x    = sosgetsol(prog, H_x);
barx22 = double(sosgetsol(prog, barx22));

fprintf('\n--- Dissipativity synthesis: certificates reported in the paper ---\n');
fprintf('P_i =\n');      disp(num2str(P, '  %.4f'));
fprintf('barZ22_i =\n'); disp(num2str(barx22, '  %.4f'));
fprintf('(eig(P) = [%s], eig(barx22) = [%s])\n\n', ...
        num2str(eig(P).', '%.4g  '), num2str(eig(barx22).', '%.4g  '));

x11 = barx11;
x12 = barx21.' * P;
x21 = x12.';
x22 = P * barx22 * P;

lmi_maxeig = NaN(size(Nvals));
for idx = 1:numel(Nvals)
    N   = Nvals(idx);
    S   = triu(ones(N), 1) - tril(ones(N), -1);
    Mn  = kron(S, eye(n));
    Z11 = kron(eye(N), x11);
    Z12 = kron(eye(N), x12);
    Z22 = kron(eye(N), x22);
    Q   = Mn.'*Z11*Mn + Mn.'*Z12 + Z12.'*Mn + Z22;
    Q   = (Q + Q.')/2;
    if N <= 200
        lmi_maxeig(idx) = max(eig(Q));
    else
        lmi_maxeig(idx) = eigs(Q, 1, 'largestreal');
    end
end

sg_margin = NaN(size(Nvals));
for idx = 1:numel(Nvals)
    N = Nvals(idx);
    arrayD = kappa * repmat(eye(3), 1, N-1);

    rng(1)
    initial = [2; 2; 2]; states = initial.'; statesdot = []; N0Tsg = [];
    usg = []; Wsg = [];
    for i = 1:T
        x0  = initial;
        M   = [x0(1), x0(2), x0(3), x0(1)*x0(3), x0(1)*x0(2)];
        M_c = [x0(1), x0(2), x0(3), x0(1)*x0(3), x0(1)*x0(2), ...
               x0(2)*x0(3), x0(2)^2, x0(1)^2, x0(3)^2];
        usg(i, :) = 100 * (-1 + 2*rand(1, 3));
        wrow      = 10/(N-1) * (-1 + 2*rand(1, 3*(N-1)));
        Wsg(i, :) = wrow;
        odeSystem = @(t, xx) A*M.' + B*usg(i, :).' + arrayD*wrow.';
        [~, X]    = ode45(odeSystem, [(i-1)*tau, i*tau], x0);
        initial   = X(end, :).';
        states    = [states; initial.'];
        statesdot = [statesdot; (A*M.' + B*usg(i, :).' + arrayD*wrow.').'];
        N0Tsg     = [N0Tsg; M_c];
    end
    U0Tsg = usg.'; W0Tsg = Wsg.';
    X1Tsg = statesdot.' + (-0.1 + 0.2*rand(3, T));
    N0Tsg = N0Tsg.';

    clear x prog
    pvar x_1 x_2 x_3
    x = [x_1; x_2; x_3];
    prog = sosprogram(x);
    b_x = [(x_1 - X_bound(1, 1)) * (X_bound(1, 2) - x_1);
           (x_2 - X_bound(2, 1)) * (X_bound(2, 2) - x_2);
           (x_3 - X_bound(3, 1)) * (X_bound(3, 2) - x_3)];
    [prog, l1] = sospolyvar(prog, monomials(x, 0:2), 'wscoeff');
    [prog, l2] = sospolyvar(prog, monomials(x, 0:2), 'wscoeff');
    [prog, l3] = sospolyvar(prog, monomials(x, 0:2), 'wscoeff');
    lam = [l1; l2; l3];
    [prog, PIsg]  = sospolyvar(prog, monomials(x, 0));
    [prog, musg]  = sospolyvar(prog, monomials(x, 0));
    [prog, invPs] = sospolymatrixvar(prog, monomials(x, 0), [n n], 'symmetric');
    [prog, Hs]    = sospolymatrixvar(prog, monomials(x, 0:1), [T n]);

    Gs = - epsilon * invPs ...
         - ((X1Tsg - arrayD*W0Tsg) * Hs) - (Hs.' * (X1Tsg - arrayD*W0Tsg).') ...
         - musg * (T * bar_varkappa * eye(n)) - PIsg * eye(n);
    Cond = [Gs, Hs.'; Hs, musg*eye(T)] - lam.' * b_x * eye(T + n);

    for ii = 1:length(lam), prog = sosineq(prog, lam(ii)); end
    prog = soseq(prog, N0Tsg * Hs - Upsilon_x * invPs);
    prog = sosmatrixineq(prog, invPs - 1e-6 * eye(n), 'quadraticMineq');
    prog = sosmatrixineq(prog, Cond, 'quadraticMineq');
    prog = sosineq(prog, PIsg - 1e-6);
    prog = sosineq(prog, musg - 1e-6);
    prog = sossolve(prog, ops);
    PIsg = double(sosgetsol(prog, PIsg));
    Ps   = double(sosgetsol(prog, invPs))^(-1);

    clear prog
    prog = sosprogram(x);
    [prog, lp1] = sospolyvar(prog, monomials(x, 0:2), 'wscoeff');
    [prog, lp2] = sospolyvar(prog, monomials(x, 0:2), 'wscoeff');
    [prog, lp3] = sospolyvar(prog, monomials(x, 0:2), 'wscoeff');
    lamp = [lp1; lp2; lp3];
    [prog, phis] = sospolyvar(prog, monomials(x, 0));
    Cond1 = (x.' * Ps * x) - phis*(x_1^2 + x_2^2 + x_3^2) - lamp.' * b_x;
    for ii = 1:length(lamp), prog = sosineq(prog, lamp(ii)); end
    prog = sosineq(prog, Cond1);
    prog = sosineq(prog, phis);
    prog = sossolve(prog, ops);
    Phis = double(sosgetsol(prog, phis));
    fprintf('N = %4d: PI_sg = %.4g, phi = %.4g, PI*phi = %.4g\n', ...
            N, PIsg, Phis, PIsg*Phis);

    Delta = (1/PIsg) * (1/Phis) * norm(arrayD)^2;
    Comp  = (ones(N) - eye(N)) * Delta;
    Composition    = ones(1, N) * (Comp - diag(repmat(epsilon, 1, N)));
    sg_margin(idx) = max(Composition);
end

fprintf('\nkappa = %g, epsilon = %g\n', kappa, epsilon);
fprintf('%8s | %32s | %32s\n', 'N', 'dissipativity LMI: max eig', ...
        'small-gain (23b): max margin');
fprintf('%8s | %32s | %32s\n', '', '(<= 0 means the LMI holds)', ...
        '(< 0 means (23b) holds)');
for idx = 1:numel(Nvals)
    fprintf('%8d | %32.4e | %32.4e\n', Nvals(idx), lmi_maxeig(idx), ...
            sg_margin(idx));
end
toc
