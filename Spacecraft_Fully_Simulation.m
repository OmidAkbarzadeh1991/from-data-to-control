clc
clear
close all
rng(1)
control_on = 1;
j1 = 20; j2 = 20; j3 = 30;
num_states = 3;
a  = [(j2 - j3)/j1; (j3 - j1)/j2; (j1 - j2)/j3];
Jv = [j1; j2; j3];
coupling_gain = 0.02;
N_sim = 2000;
num_runs = 100;
tau   = 1e-3;
tspan = 20;
T     = tspan / tau;
States = zeros(T + 1, num_states, num_runs);
for r = 1:num_runs
    Xmat = -2 + 4 * rand(num_states, N_sim);
    States(1, :, r) = Xmat(:, 1).';
    for i = 1:T
        X1 = Xmat(1, :);  X2 = Xmat(2, :);  X3 = Xmat(3, :);
        u1 = -92.3296*X1.^2 - 42.4488*X1.*X2 + 145.568*X1.*X3 + 84.2951*X2.^2 ...
             - 489.389*X2.*X3 + 156.6679*X3.^2 - 38251.1977*X1 - 315.5688*X2 ...
             + 2551.2783*X3;
        u2 = 143.3514*X1.^2 - 776.2086*X1.*X2 + 360.0186*X1.*X3 + 251.1473*X2.^2 ...
             + 722.9704*X2.*X3 - 160.4184*X3.^2 + 212.894*X1 - 15720.342*X2 ...
             - 10594.1732*X3;
        u3 = -103.8397*X1.^2 + 1255.1938*X1.*X2 - 647.8731*X1.*X3 - 554.7252*X2.^2 ...
             - 316.4707*X2.*X3 + 94.4806*X3.^2 + 1531.45*X1 - 15774.9915*X2 ...
             - 44634.7003*X3;
        f = [a(1) * X2.*X3;
             a(2) * X3.*X1;
             a(3) * X1.*X2] + control_on * [u1; u2; u3] ./ Jv;
        f = f + coupling_gain * (sum(Xmat, 2) - Xmat) ./ Jv;
        Xmat = Xmat + tau * f;
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
S_i = @(x1, x2, x3) 8.2063.*x1.^2 + 0.09464.*x1.*x2 - 1.0109.*x1.*x3 ...
                  + 3.4732.*x2.^2 + 4.3726.*x2.*x3 + 6.4915.*x3.^2;
eta_i = 108.3956;
mu_i  = 110.8317;
x_range = linspace(-10, 10, 400);
[X1, X2, X3] = meshgrid(x_range, x_range, x_range);
V_eta = S_i(X1, X2, X3) - eta_i;
V_mu  = S_i(X1, X2, X3) - mu_i;
s1 = isosurface(X1, X2, X3, V_eta, 0);
p1 = patch(s1);
reducepatch(p1, 0.2);
isonormals(X1, X2, X3, V_eta, p1);
set(p1, 'FaceColor', 'green', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
s2 = isosurface(X1, X2, X3, V_mu, 0);
p2 = patch(s2);
reducepatch(p2, 0.2);
isonormals(X1, X2, X3, V_mu, p2);
set(p2, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
X0_bound  = 2 * [-1 1; -1 1; -1 1];
X1_bound1 = [ 2.5 5; -5 -3; -5 -4];
X1_bound2 = [ 3.5 5;  3.5 5; -2 2];
X1_bound3 = [-5 -4;   4  5; 2.5 5];
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
    plotFace(verts0, faces(i, :), 'blue');
    plotFace(verts1, faces(i, :), 'red');
    plotFace(verts2, faces(i, :), 'red');
    plotFace(verts3, faces(i, :), 'red');
end
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
