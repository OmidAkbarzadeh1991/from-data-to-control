clc
clear
close all

ops.solver = 'mosek';

rng(1)

tic

X_bound   = [-5 5; -5 5; -5 5];
X0_bound  = 2 * [-1 1; -1 1; -1 1];
X1_bound1 = [ 2.5 5; -5 -3; -5 -4];
X1_bound2 = [ 3.5 5;  3.5 5; -2 2];
X1_bound3 = [-5 -4;   4  5; 2.5 5];

number_subsystems = 2000;

T   = 17;
tau = 0.05;

j1 = 20;
j2 = 20;
j3 = 30;

initial = -0.1 * [0.5; 0.5; 0.5];

A = blkdiag((j2 - j3)/j1, (j3 - j1)/j2, (j1 - j2)/j3);
A = [zeros(3, 3), A];

B = blkdiag(1/j1, 1/j2, 1/j3);

D_paper = 0.02 * B;
kappa   = 0.02;
D       = D_paper / kappa;

Matrix_D = D;

states    = initial.';
statesdot = [];
N0T       = [];
Vector_W  = [];

for i = 1:T

    x0 = initial;

    M  = [x0(1), x0(2), x0(3), x0(2)*x0(3), x0(3)*x0(1), x0(1)*x0(2)];

    M1 = [x0(1), x0(2), x0(3), x0(2)*x0(3), x0(3)*x0(1), x0(1)*x0(2), ...
          x0(1)^2, x0(2)^2, x0(3)^2];

    u(:, i) = 3 * (-100 + 200 * rand(3, 1));

    w = 0.08 * [-10 + 20*rand, -10 + 20*rand, -20 + 20*rand];
    Vector_W = vertcat(Vector_W, w);

    odeSystem = @(t, x) A*M.' + B*u(:, i) + Matrix_D * Vector_W(i, :).';
    tspan     = [(i-1)*tau, i*tau];
    options   = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
    [~, X]    = ode45(odeSystem, tspan, x0, options);

    initial   = X(end, :).';
    states    = [states; initial.'];
    statesdot = [statesdot; (A*M.' + B*u(:, i) + Matrix_D * Vector_W(i, :).').'];
    N0T       = [N0T; M1];

end

U0T = u;
W0T = Vector_W.';
X0T = states(1:end-1, :).';

X1T = statesdot.' + (-0.5 + rand(3, T));

N0T = N0T.';

clear tau initial states statesdot x0 u X i w odeSystem tspan options

pvar x_1 x_2 x_3
x = [x_1; x_2; x_3];

Upsilon_x = [1    0    0;
             0    1    0;
             0    0    1;
             0    0    x_2;
             x_3  0    0;
             x_2  0    0;
             x_1  0    0;
             0    x_2  0;
             0    0    x_3];

prog = sosprogram(x);

b_x_1 = (x_1 - X_bound(1, 1)) * (X_bound(1, 2) - x_1);
b_x_2 = (x_2 - X_bound(2, 1)) * (X_bound(2, 2) - x_2);
b_x_3 = (x_3 - X_bound(3, 1)) * (X_bound(3, 2) - x_3);
b_x   = [b_x_1; b_x_2; b_x_3];

Mon_Lagrangian = monomials(x, 0:2);

[prog, L3_1] = sospolyvar(prog, Mon_Lagrangian, 'wscoeff');
[prog, L3_2] = sospolyvar(prog, Mon_Lagrangian, 'wscoeff');
[prog, L3_3] = sospolyvar(prog, Mon_Lagrangian, 'wscoeff');
lambda_x = [L3_1; L3_2; L3_3];

PI           = 0.6356;
epsilon      = 0.99;
bar_varkappa = 3 * 0.25;

[prog, mu]   = sospolyvar(prog, monomials(x, 0));
[prog, invP] = sospolymatrixvar(prog, monomials(x, 0), ...
                                [length(x), length(x)], 'symmetric');
[prog, H_x]  = sospolymatrixvar(prog, monomials(x, 0:1), [T, length(x)]);

[prog, barx22] = sospolymatrixvar(prog, monomials(x, 0), ...
                                  [length(x), length(x)], 'symmetric');

barx11 = (1/PI) * (Matrix_D.' * Matrix_D);
barx21 = zeros(length(x));

G = - epsilon * invP ...
    - ((X1T - Matrix_D * W0T) * H_x) - (H_x.' * (X1T - Matrix_D * W0T).') ...
    - mu * (T * bar_varkappa * eye(length(x))) ...
    - PI * eye(length(x));

Condition_3_creation = ...
    [ G + barx22,  H_x.';
      H_x,         mu * eye(T) ] ...
    - lambda_x.' * b_x * eye(T + length(x));

for ii = 1:length(lambda_x)
    prog = sosineq(prog, lambda_x(ii));
end

prog = soseq(prog, N0T * H_x - Upsilon_x * invP);

prog = sosmatrixineq(prog, invP - 1e-6 * eye(length(x)), 'quadraticMineq');

prog = sosmatrixineq(prog, Condition_3_creation, 'quadraticMineq');

prog = sosineq(prog, mu - 1e-6);

prog = sosmatrixineq(prog, -barx22, 'quadraticMineq');

prog = sossolve(prog, ops);

mu   = double(sosgetsol(prog, mu));
invP = double(sosgetsol(prog, invP));
P    = invP^(-1);
H_x  = sosgetsol(prog, H_x);

barx22 = double(sosgetsol(prog, barx22));

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
disp('barx22 ='); disp(barx22);

RT = toc;
fprintf('Per-subsystem runtime RT = %.2f seconds\n', RT);

allvars = whos;
MU_bytes = sum([allvars.bytes]);
fprintf('Per-subsystem memory  MU = %.2f MB\n', MU_bytes / 1e6);

n_i = 3;
Nn  = number_subsystems * n_i;

x11 = barx11;
x12 = barx21.' * P;
x21 = x12.';
x22 = P * barx22 * P;

connection = kron(ones(number_subsystems) - eye(number_subsystems), kappa * eye(n_i));

X1A = kron(eye(number_subsystems), x11);
X2A = kron(eye(number_subsystems), x12);
X3A = X2A.';
X4A = kron(eye(number_subsystems), x22);

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
eig_lmi = real(eig((Dissipativity + Dissipativity.')/2));

if all(eig_lmi < 0)
    fprintf('min eig = %.4e,  max eig = %.4e\n', min(eig_lmi), max(eig_lmi));
    msg = 'Compositionality condition satisfied';
else
    fprintf('max eig = %.4e (should be negative)\n', max(eig_lmi));
    msg = 'Compositionality condition NOT satisfied';
end
border = repmat('-', 1, length(msg) + 4);
disp(border); disp(['* ', msg, ' *']); disp(border);

clear x prog X1A X2A X3A X4A Comp connection Dissipativity

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

gamma_network = number_subsystems * gamma
beta_network  = number_subsystems * beta
