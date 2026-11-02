clc
clear
close all

ops.solver = 'mosek';

rng(1)

tic

X_bound   = [-20 20; -20 20; -20 20];
X0_bound  = 3 * [-1 1; -1 1; -1 1];
X1_bound1 = [-20 -4; -20 -15;  4  20];
X1_bound2 = [  5 20;  11  20;  4  20];
X1_bound3 = [  5 20;  11  20; -20 -4];

Subsystems = 1000;

T   = 15;
tau = 0.01;

sigma = 10;
r     = 28;
b     = 8/3;

initial = [2; 2; 2];

A = [-sigma, sigma,  0,  0, 0;
      r,    -1,      0, -1, 0;
      0,     0,     -b,  0, 1];

B = eye(3);

D_coeff = 1;
D = -D_coeff * eye(3);

states    = initial.';
statesdot = [];
N0T       = [];

for i = 1:T

    x0 = initial;

    M   = [x0(1), x0(2), x0(3), x0(1)*x0(3), x0(1)*x0(2)];

    M_c = [x0(1), x0(2), x0(3), x0(1)*x0(3), x0(1)*x0(2), x0(2)*x0(3), ...
           x0(2)^2, x0(1)^2, x0(3)^2];

    u1(i, :) = 100 * (-1 + 2 * rand);
    u2(i, :) = 100 * (-1 + 2 * rand);
    u3(i, :) = 100 * (-1 + 2 * rand);
    u(i, :)  = [u1(i, :), u2(i, :), u3(i, :)];

    w(i, :) = 999 * 0.001 * [-10 + 20*rand, -10 + 20*rand, -10 + 20*rand];

    odeSystem = @(t, x) A * M.' + B * u(i, :).' + D * w(i, :).';
    tspan     = [(i-1)*tau, i*tau];
    [~, X]    = ode45(odeSystem, tspan, x0);

    initial   = X(end, :).';
    states    = [states; initial.'];
    statesdot = [statesdot; (A*M.' + B*u(i, :).' + D*w(i, :).').'];
    N0T       = [N0T; M_c];

end

U0T = u.';
W0T = w.';
X0T = states(1:end-1, :).';

X1T = statesdot.' + (-0.1 + 0.2 * rand(3, T));

N0T = N0T.';

clear tau initial states statesdot x0 u u1 u2 u3 w X i odeSystem tspan

pvar x_1 x_2 x_3
x = [x_1; x_2; x_3];

Upsilon_x = [1    0    0;
             0    1    0;
             0    0    1;
             x_3  0    0;
             x_2  0    0;
             0    x_3  0;
             0    x_2  0;
             x_1  0    0;
             0    0    x_3];

prog = sosprogram(x);

b_x_1 = (x_1 - X_bound(1, 1)) * (X_bound(1, 2) - x_1);
b_x_2 = (x_2 - X_bound(2, 1)) * (X_bound(2, 2) - x_2);
b_x_3 = (x_3 - X_bound(3, 1)) * (X_bound(3, 2) - x_3);
b_x   = [b_x_1; b_x_2; b_x_3];

Mon_Lagrangian = monomials(x, 0:2);

[prog, lambda_x_1] = sospolyvar(prog, Mon_Lagrangian, 'wscoeff');
[prog, lambda_x_2] = sospolyvar(prog, Mon_Lagrangian, 'wscoeff');
[prog, lambda_x_3] = sospolyvar(prog, Mon_Lagrangian, 'wscoeff');
lambda_x = [lambda_x_1; lambda_x_2; lambda_x_3];

PI           = 0.8481;
epsilon      = 0.99;
bar_varkappa = 3 * 0.01;

[prog, mu]   = sospolyvar(prog, monomials(x, 0));
[prog, invP] = sospolymatrixvar(prog, monomials(x, 0), ...
                                [length(x), length(x)], 'symmetric');
[prog, H_x]  = sospolymatrixvar(prog, monomials(x, 0:1), [T, length(x)]);

[prog, barx11] = sospolymatrixvar(prog, monomials(x, 0), ...
                                  [length(x), length(x)], 'symmetric');
[prog, barx22] = sospolymatrixvar(prog, monomials(x, 0), ...
                                  [length(x), length(x)], 'symmetric');
[prog, barx21] = sospolymatrixvar(prog, monomials(x, 0), ...
                                  [length(x), length(x)], 'symmetric');

G = - epsilon * invP ...
    - ((X1T - D * W0T) * H_x) - (H_x.' * (X1T - D * W0T).') ...
    - mu * (T * bar_varkappa * eye(length(x))) ...
    - PI * eye(length(x));

Condition_3_creation = ...
    [ G + barx22,  H_x.',          barx21;
      H_x,         mu * eye(T),    zeros(T, length(x));
      barx21.',    zeros(T, length(x)).', barx11 - (1/PI) * D.' * D ] ...
    - lambda_x.' * b_x * eye(T + 2 * length(x));

for ii = 1:length(lambda_x)
    prog = sosineq(prog, lambda_x(ii));
end

prog = soseq(prog, N0T * H_x - Upsilon_x * invP);

prog = sosmatrixineq(prog, invP - 1e-6 * eye(length(x)), 'quadraticMineq');

prog = sosmatrixineq(prog, Condition_3_creation, 'quadraticMineq');

prog = sosineq(prog, mu - 1e-6);

prog = sosineq(prog,  barx11);
prog = soseq (prog,  barx11 - barx11.');
prog = soseq (prog,  barx22 - barx22.');
prog = sosineq(prog, -barx22);
prog = sosineq(prog, -barx21);

prog = sossolve(prog, ops);

mu   = double(sosgetsol(prog, mu));
invP = double(sosgetsol(prog, invP));
P    = invP^(-1);
H_x  = sosgetsol(prog, H_x);

barx11 = double(sosgetsol(prog, barx11));
barx22 = double(sosgetsol(prog, barx22));
barx21 = double(sosgetsol(prog, barx21));

msg = 'Control input & storage certificate (CSC)';
border = repmat('-', 1, length(msg) + 4);
disp(border); disp(['* ', msg, ' *']); disp(border);

Control_input       = U0T * H_x * P * x
Storage_certificate = x.' * P * x

fprintf('pi_i  = %.4f\n', PI);
fprintf('mu    = %.4f (noise multiplier)\n', mu);
disp('P_i ='); disp(P);
disp('Eigenvalues of P_i:'); disp(eig(P));

eig_barx22 = eig((barx22 + barx22.')/2);
fprintf('max eig(barx22) = %.3e  ->  barx22 <= 0: %s\n', max(eig_barx22), ...
        string(max(eig_barx22) <= 1e-8));
eig_x11gap = eig(barx11 - (1/PI) * (D.' * D));
fprintf('min eig(barx11 - (1/pi)D''D) = %.3e  ->  satisfied: %s\n', ...
        min(eig_x11gap), string(min(eig_x11gap) >= -1e-8));

RT = toc;
fprintf('Per-subsystem runtime RT = %.2f seconds\n', RT);

allvars = whos;
MU_bytes = sum([allvars.bytes]);
fprintf('Per-subsystem memory  MU = %.2f MB\n', MU_bytes / 1e6);

n_i = 3;
Nn  = Subsystems * n_i;

connection = kron(ones(Subsystems) - eye(Subsystems), 1e-5 * eye(n_i));

x11 = barx11;
x12 = barx21.' * P;
x21 = x12.';
x22 = P * barx22 * P;

X1A = kron(eye(Subsystems), x11);
X2A = kron(eye(Subsystems), x12);
X3A = X2A.';
X4A = kron(eye(Subsystems), x22);

Comp = [X1A, X2A; X3A, X4A];

disp('Checking the symmetry of the supply-rate matrix:');
if issymmetric(round(Comp, 10))
    msg = 'Supply-rate matrix is symmetric';
else
    msg = 'Supply-rate matrix is NOT symmetric';
end
border = repmat('-', 1, length(msg) + 4);
disp(border); disp(['* ', msg, ' *']); disp(border);

Dissipativity = [connection; eye(Nn)].' * Comp * [connection; eye(Nn)];

disp('Checking the dissipativity (compositionality) condition:');
eig_lmi = real(eig(Dissipativity));

if all(eig_lmi < 0)
    fprintf('min eig = %.4e,  max eig = %.4e\n', min(eig_lmi), max(eig_lmi));
    msg = 'Compositionality condition satisfied';
else
    fprintf('max eig = %.4e (should be negative)\n', max(eig_lmi));
    msg = 'Compositionality condition NOT satisfied';
end
border = repmat('-', 1, length(msg) + 4);
disp(border); disp(['* ', msg, ' *']); disp(border);

clear x prog X1A X2A X3A X4A

pvar x1 x2 x3
var_x = [x1; x2; x3];

prog = sosprogram(var_x);

b0_1 = (x1 - X0_bound(1, 1)) * (X0_bound(1, 2) - x1);
b0_2 = (x2 - X0_bound(2, 1)) * (X0_bound(2, 2) - x2);
b0_3 = (x3 - X0_bound(3, 1)) * (X0_bound(3, 2) - x3);
b0   = [b0_1; b0_2; b0_3];

b_u_1 = (x1 - X1_bound1(1, 1)) * (X1_bound1(1, 2) - x1);
b_u_2 = (x2 - X1_bound1(2, 1)) * (X1_bound1(2, 2) - x2);
b_u_3 = (x3 - X1_bound1(3, 1)) * (X1_bound1(3, 2) - x3);
b_u   = [b_u_1; b_u_2; b_u_3];

b_uu_1 = (x1 - X1_bound2(1, 1)) * (X1_bound2(1, 2) - x1);
b_uu_2 = (x2 - X1_bound2(2, 1)) * (X1_bound2(2, 2) - x2);
b_uu_3 = (x3 - X1_bound2(3, 1)) * (X1_bound2(3, 2) - x3);
b_uu   = [b_uu_1; b_uu_2; b_uu_3];

b_uuu_1 = (x1 - X1_bound3(1, 1)) * (X1_bound3(1, 2) - x1);
b_uuu_2 = (x2 - X1_bound3(2, 1)) * (X1_bound3(2, 2) - x2);
b_uuu_3 = (x3 - X1_bound3(3, 1)) * (X1_bound3(3, 2) - x3);
b_uuu   = [b_uuu_1; b_uuu_2; b_uuu_3];

Mon_Lagrangian1 = monomials(var_x, 0:2);

[prog, L0_1] = sospolyvar(prog, Mon_Lagrangian1, 'wscoeff');
[prog, L0_2] = sospolyvar(prog, Mon_Lagrangian1, 'wscoeff');
[prog, L0_3] = sospolyvar(prog, Mon_Lagrangian1, 'wscoeff');
Lambda0 = [L0_1; L0_2; L0_3];

[prog, L1_1] = sospolyvar(prog, Mon_Lagrangian1, 'wscoeff');
[prog, L1_2] = sospolyvar(prog, Mon_Lagrangian1, 'wscoeff');
[prog, L1_3] = sospolyvar(prog, Mon_Lagrangian1, 'wscoeff');
Lambda_u_1 = [L1_1; L1_2; L1_3];

[prog, L2_1] = sospolyvar(prog, Mon_Lagrangian1, 'wscoeff');
[prog, L2_2] = sospolyvar(prog, Mon_Lagrangian1, 'wscoeff');
[prog, L2_3] = sospolyvar(prog, Mon_Lagrangian1, 'wscoeff');
Lambda_u_2 = [L2_1; L2_2; L2_3];

[prog, L4_1] = sospolyvar(prog, Mon_Lagrangian1, 'wscoeff');
[prog, L4_2] = sospolyvar(prog, Mon_Lagrangian1, 'wscoeff');
[prog, L4_3] = sospolyvar(prog, Mon_Lagrangian1, 'wscoeff');
Lambda_u_3 = [L4_1; L4_2; L4_3];

Barrier = var_x.' * P * var_x;

[prog, gamma] = sospolyvar(prog, monomials(var_x, 0));
[prog, beta]  = sospolyvar(prog, monomials(var_x, 0));

prog = sosineq(prog, beta - gamma);
prog = sosineq(prog, gamma);
prog = sosineq(prog, beta);

prog = sosineq(prog, -Barrier - Lambda0.'    * b0    + gamma);
prog = sosineq(prog,  Barrier - Lambda_u_1.' * b_u   - beta);
prog = sosineq(prog,  Barrier - Lambda_u_2.' * b_uu  - beta);
prog = sosineq(prog,  Barrier - Lambda_u_3.' * b_uuu - beta);

for jj = 1:length(Lambda0),    prog = sosineq(prog, Lambda0(jj));    end
for kk = 1:length(Lambda_u_1), prog = sosineq(prog, Lambda_u_1(kk)); end
for zz = 1:length(Lambda_u_2), prog = sosineq(prog, Lambda_u_2(zz)); end
for cc = 1:length(Lambda_u_3), prog = sosineq(prog, Lambda_u_3(cc)); end

prog = sossolve(prog, ops);

SOLV0 = sosgetsol(prog, -Barrier - Lambda0.'    * b0    + gamma);
SOLV1 = sosgetsol(prog,  Barrier - Lambda_u_1.' * b_u   - beta);
SOLV2 = sosgetsol(prog,  Barrier - Lambda_u_2.' * b_uu  - beta);
SOLV3 = sosgetsol(prog,  Barrier - Lambda_u_3.' * b_uuu - beta);

[Q0, Z0] = findsos(SOLV0, ops);
[Q1, Z1] = findsos(SOLV1, ops);
[Q2, Z2] = findsos(SOLV2, ops);
[Q3, Z3] = findsos(SOLV3, ops);

sos_ok = ~[isempty(Q0), isempty(Q1), isempty(Q2), isempty(Q3)];
fprintf('Independent SOS verification (initial, unsafe 1-3): %d %d %d %d\n', sos_ok);

msg = 'Level sets';
border = repmat('-', 1, length(msg) + 4);
disp(border); disp(['* ', msg, ' *']); disp(border);

gamma = double(sosgetsol(prog, gamma))
beta  = double(sosgetsol(prog, beta))

fprintf('Strict separation mu_i > eta_i: %s (gap = %.4f)\n', ...
        string(beta > gamma), beta - gamma);

gamma_network = Subsystems * gamma
beta_network  = Subsystems * beta
