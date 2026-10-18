%% ========================================================================
%  Closed-Loop Simulation and Plots
%  Benchmark: Spacecraft Network with Star Interconnection Topology
%
%  Companion code for:
%    "From Data to Control: A Formal Compositional Framework for
%     Large-Scale Interconnected Networks" (arXiv:2409.12469)
%
%  Simulates the network under the data-driven decentralized controller
%  synthesized by Spacecraft_Star_Network.m and plots closed-loop
%  trajectories with the certificate level sets and the regions of interest.
%
%  In a star topology each leaf depends only on its own state and the hub,
%  so simulating N_sim subsystems (hub + N_sim-1 leaves) reproduces EXACTLY
%  the trajectories those subsystems have inside the full 2000-node network.
%% ========================================================================

clc
clear
close all

rng(1)                          % fix random initial conditions

%% ===================== Setup ============================================

control_on = 1;                 % 1: closed loop, 0: open loop

j1 = 20; j2 = 20; j3 = 30;      % principal moments of inertia

num_states = 3;

% subsystem matrices acting on R(x) = [x1; x2; x3; x2x3; x3x1; x1x2]
A = [0 0 0 (j2 - j3)/j1 0 0;
     0 0 0 0 (j3 - j1)/j2 0;
     0 0 0 0 0 (j1 - j2)/j3];

B = blkdiag(1/j1, 1/j2, 1/j3);

coupling_gain = 0.02;           % star coupling D = 0.02*B (as in the paper)
D = coupling_gain * B;

N_sim = 2000;                    % simulated subsystems (exact subset of the
                                % 2000-node star; see header)

x = -2 + 4 * rand(num_states * N_sim, 1);   % initial states in X0 = [-2,2]^3

%% ===================== Network assembly =================================

m = 3;                          % rows of A
n = 6;                          % columns of A

A_network = zeros(N_sim*m, N_sim*n);
B_network = kron(eye(N_sim), B);

for i = 1:N_sim
    A_network((i-1)*m + 1 : i*m, (i-1)*n + 1 : i*n) = A;
    if i < N_sim
        A_network(i*m + 1 : (i+1)*m, 1:3) = D;   % leaf i+1 <- hub states
    end
end

%% ===================== Controller and monomials =========================
%  Coefficients as synthesized by Spacecraft_Star_Network.m (paper values).

syms x_1 x_2 x_3

u_symbolic = control_on * [
    -92.3296*x_1^2 - 42.4488*x_1*x_2 + 145.568*x_1*x_3 + 84.2951*x_2^2 ...
    - 489.389*x_2*x_3 + 156.6679*x_3^2 - 38251.1977*x_1 - 315.5688*x_2 ...
    + 2551.2783*x_3;
     143.3514*x_1^2 - 776.2086*x_1*x_2 + 360.0186*x_1*x_3 + 251.1473*x_2^2 ...
    + 722.9704*x_2*x_3 - 160.4184*x_3^2 + 212.894*x_1 - 15720.342*x_2 ...
    - 10594.1732*x_3;
    -103.8397*x_1^2 + 1255.1938*x_1*x_2 - 647.8731*x_1*x_3 - 554.7252*x_2^2 ...
    - 316.4707*x_2*x_3 + 94.4806*x_3^2 + 1531.45*x_1 - 15774.9915*x_2 ...
    - 44634.7003*x_3];

M_x = [x_1; x_2; x_3; x_2*x_3; x_3*x_1; x_1*x_2];

% fast numerical evaluation
controller_fast = matlabFunction(u_symbolic, 'Vars', {[x_1; x_2; x_3]});
monomial_fast   = matlabFunction(M_x,        'Vars', {[x_1; x_2; x_3]});

%% ===================== Simulation (forward Euler) =======================
%  tau = 1e-3 is below the Euler stability limit of the fastest closed-loop
%  mode (~ -38251/20 per second, limit ~ 1.05e-3); reduce tau if ringing
%  appears near t = 0.

tau   = 1e-3;                   % integration step
tspan = 20;                     % simulation time
T     = tspan / tau;

States = zeros(T + 1, num_states * N_sim);
States(1, :) = x.';

for i = 1:T

    u    = zeros(num_states * N_sim, 1);
    mono = zeros(n * N_sim, 1);

    for j = 1:N_sim
        idx     = (j-1)*3 + 1;
        x_local = x(idx : idx+2);

        u((j-1)*3 + 1 : j*3)   = controller_fast(x_local);
        mono((j-1)*n + 1 : j*n) = monomial_fast(x_local);
    end

    x = x + tau * (A_network * mono + B_network * u);

    States(i + 1, :) = x.';

end

%% ===================== Plots ============================================

num_plots = 120;                % subsystems to plot

figure
hold on;

random_subsystems = randperm(N_sim, num_plots);

for i = 1:num_plots
    sub_idx = random_subsystems(i);
    plot3(States(:, 3*sub_idx - 2), States(:, 3*sub_idx - 1), ...
          States(:, 3*sub_idx), 'LineWidth', 2, 'LineStyle', '-');
end

% ----- per-subsystem certificate level sets {S_i = eta_i}, {S_i = mu_i} ----
S_i = @(x1, x2, x3) 8.2063.*x1.^2 + 0.09464.*x1.*x2 - 1.0109.*x1.*x3 ...
                  + 3.4732.*x2.^2 + 4.3726.*x2.*x3 + 6.4915.*x3.^2;

eta_i = 127.7789;               % initial level set
mu_i  = 132.5868;               % unsafe level set

x_range = linspace(-10, 10, 200);
[X1, X2, X3] = meshgrid(x_range, x_range, x_range);

V_eta = S_i(X1, X2, X3) - eta_i;
V_mu  = S_i(X1, X2, X3) - mu_i;

s1 = isosurface(X1, X2, X3, V_eta, 0);        % {S_i = eta_i} in green
p1 = patch(s1);
isonormals(X1, X2, X3, V_eta, p1);
set(p1, 'FaceColor', 'green', 'EdgeColor', 'none', 'FaceAlpha', 0.3);

s2 = isosurface(X1, X2, X3, V_mu, 0);         % {S_i = mu_i} in red
p2 = patch(s2);
isonormals(X1, X2, X3, V_mu, p2);
set(p2, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% ----- regions of interest -------------------------------------------------
X0_bound  = 2 * [-1 1; -1 1; -1 1];    % initial set
X1_bound1 = [ 2.5 5; -5 -3; -5 -4];    % unsafe box 1
X1_bound2 = [ 2.5 5;  4  5; 2.5 5];    % unsafe box 2
X1_bound3 = [-5 -4;   4  5; 2.5 5];    % unsafe box 3

calcVertices = @(b) [b(1,1), b(2,1), b(3,1); ...
                     b(1,2), b(2,1), b(3,1); ...
                     b(1,2), b(2,2), b(3,1); ...
                     b(1,1), b(2,2), b(3,1); ...
                     b(1,1), b(2,1), b(3,2); ...
                     b(1,2), b(2,1), b(3,2); ...
                     b(1,2), b(2,2), b(3,2); ...
                     b(1,1), b(2,2), b(3,2)];

plotFace = @(v, f, c) fill3(v(f, 1), v(f, 2), v(f, 3), c, ...
                            'FaceAlpha', 0.4, 'EdgeColor', 'none');

verts0 = calcVertices(X0_bound);
verts1 = calcVertices(X1_bound1);
verts2 = calcVertices(X1_bound2);
verts3 = calcVertices(X1_bound3);

faces = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];

for i = 1:size(faces, 1)
    plotFace(verts0, faces(i, :), 'blue');    % initial set (blue)
    plotFace(verts1, faces(i, :), 'red');     % unsafe box 1 (red)
    plotFace(verts2, faces(i, :), 'red');     % unsafe box 2 (red)
    plotFace(verts3, faces(i, :), 'red');     % unsafe box 3 (red)
end

% ----- axes and styling ------------------------------------------------------
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6);
axis equal;
xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
zlabel('$x_3$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 28)
set(gca, 'GridLineStyle', ':', 'GridAlpha', 1, 'LineWidth', 1);
box on;

view(3);
xlim([-5 5]);
ylim([-5 5]);
zlim([-5 5]);