clc
clear
close all

rng(1)

control_on = 1;

sigma  = 10;
rho    = 28;
beta_L = 8/3;

N = 1000;

d = 1e-5;

num_runs = 100;

tau   = 1e-4;
tspan = 1;
T     = tspan / tau;

States = zeros(T + 1, 3, num_runs);

for r = 1:num_runs

    Xmat = 3 * (-1 + 2 * rand(3, N));
    States(1, :, r) = Xmat(:, 1).';

    for i = 1:T

        X1 = Xmat(1, :);  X2 = Xmat(2, :);  X3 = Xmat(3, :);

        u1 = -4.3235*X1.^2 + 12.7482*X1.*X2 - 3.3942*X1.*X3 - 3.5729*X2.^2 ...
             + 1.9325*X2.*X3 + 4.4843*X3.^2 - 552.0742*X1 + 33.4933*X2 ...
             + 138.1536*X3;
        u2 = -38.2107*X1.^2 + 15.2541*X1.*X2 - 3.8777*X1.*X3 - 0.68325*X2.^2 ...
             + 5.0418*X2.*X3 - 13.4643*X3.^2 + 160.3163*X1 - 165.9182*X2 ...
             + 145.0549*X3;
        u3 = 9.5672*X1.^2 + 2.8845*X1.*X2 - 19.0905*X1.*X3 - 6.0917*X2.^2 ...
             + 13.9545*X2.*X3 + 2.4699*X3.^2 - 254.5404*X1 - 140.8946*X2 ...
             - 181.1664*X3;

        f = [ -sigma*X1 + sigma*X2;
               rho*X1 - X2 - X1.*X3;
              -beta_L*X3 + X1.*X2 ];

        coupling = -d * (sum(Xmat, 2) - Xmat);

        Xmat = Xmat + tau * (f + control_on * [u1; u2; u3] + coupling);

        States(i + 1, :, r) = Xmat(:, 1).';

    end

end

figure
hold on;

k = 1:50:size(States, 1);
for r = 1:num_runs
    plot3(States(k, 1, r), States(k, 2, r), States(k, 3, r), ...
          'LineWidth', 2, 'LineStyle', '-');
end

S_i = @(x1, x2, x3) 4.6885.*x1.^2 - 1.2132.*x1.*x2 - 0.4946.*x1.*x3 ...
                  + 1.4175.*x2.^2 + 0.37465.*x2.*x3 + 1.5099.*x3.^2;

eta_i = 98.2169;
mu_i  = 100.6141;

g_range = linspace(-20, 20, 400);
[g1, g2, g3] = meshgrid(g_range, g_range, g_range);

V_eta = S_i(g1, g2, g3) - eta_i;
V_mu  = S_i(g1, g2, g3) - mu_i;

s1 = isosurface(g1, g2, g3, V_eta, 0);
p1 = patch(s1);
reducepatch(p1, 0.2);
isonormals(g1, g2, g3, V_eta, p1);
set(p1, 'FaceColor', 'green', 'EdgeColor', 'none', 'FaceAlpha', 0.3);

s2 = isosurface(g1, g2, g3, V_mu, 0);
p2 = patch(s2);
reducepatch(p2, 0.2);
isonormals(g1, g2, g3, V_mu, p2);
set(p2, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.3);

X0_bound  = 3 * [-1 1; -1 1; -1 1];
X1_bound1 = [-20 -4; -20 -15;  4  20];
X1_bound2 = [  5 20;  11  20;  4  20];
X1_bound3 = [  5 20;  11  20; -20 -4];

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
    plotFace(verts0, faces(i, :), 'b');
    plotFace(verts1, faces(i, :), 'r');
    plotFace(verts2, faces(i, :), 'r');
    plotFace(verts3, faces(i, :), 'r');
end

xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
zlabel('$x_3$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
axis equal;
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 28)
set(gca, 'GridLineStyle', ':', 'GridAlpha', 1, 'LineWidth', 1);
box on;
view(3);

xlim([-20 20]);
ylim([-20 20]);
zlim([-20 20]);
