%% ========================================================================
%  Space and Time Complexity: Compositional versus Monolithic Synthesis
%
%  Companion code for:
%    "From Data to Control: A Formal Compositional Framework for
%     Large-Scale Interconnected Networks" (arXiv:2409.12469)
%
%  Evaluates the per-subsystem complexity estimates of the paper
%  (decision-variable count kappa_i, space cost C_i, time cost J_i) and
%  compares the compositional totals N*C_i and N*J_i against the
%  monolithic synthesis, whose quantities grow with the aggregated
%  dimension N*n_i, horizon N*T, and dictionary size N*M_i.
%  Parameters follow the spacecraft benchmark. Requires only base MATLAB.
%% ========================================================================

clc
clear
close all

%% ===================== Parameters =======================================

numStatesLocal    = 3;      % n_i (subsystem dimension)
horizonLocal      = 17;     % T   (data horizon)
degreeHi          = 2;      % r_i (degree of H_i, matched to the degree-2
                            %      dictionary; the benchmark codes implement
                            %      H_i of degree 1, a uniform overestimate)
numThetaRowsLocal = 9;      % M_i (rows of Theta_i)
epsi              = 1e-8;   % solver accuracy epsilon

Nvals = 1:2000;             % number of subsystems N

%% ===================== Per-subsystem quantities =========================

% monomials in n_i variables of degree <= r_i
numMonomialsLocal = nchoosek(numStatesLocal + degreeHi, degreeHi);

% decision-variable count kappa_i:
%   T*n*binom(n+r,r)  (H_i)  +  3*n*(n+1)/2  (Zbar11, Zbar22, S_i symmetric)
%   + n^2  (Zbar12 nonsymmetric; Zbar21 = Zbar12' not free)
numDecisionVarsLocal = horizonLocal*numStatesLocal*numMonomialsLocal ...
                     + 1.5*numStatesLocal*(numStatesLocal + 1) ...
                     + numStatesLocal^2;

% interior-point iteration factor: sqrt(2n+T)*log(1/eps)
iterFactorLocal = sqrt(2*numStatesLocal + horizonLocal) * log(1/epsi);

% time cost J_i and compositional total J_comp = N*J_i
J_i_time = iterFactorLocal * (numDecisionVarsLocal^3);
J_comp   = Nvals .* J_i_time;

% space cost C_i and compositional total C_comp = N*C_i:
%   kappa*(2n+T)^2 + kappa^2 + (T*n*binom + n*(n+1)/2) * (M*n*binom)
C_i_space = numDecisionVarsLocal*(2*numStatesLocal + horizonLocal)^2 ...
          + numDecisionVarsLocal^2 ...
          + (horizonLocal*numStatesLocal*numMonomialsLocal ...
             + numStatesLocal*(numStatesLocal + 1)/2) ...
            * (numThetaRowsLocal*numStatesLocal*numMonomialsLocal);
C_comp    = Nvals .* C_i_space;

%% ===================== Monolithic quantities ============================

numStatesTotal    = Nvals * numStatesLocal;        % N*n_i
horizonTotal      = Nvals * horizonLocal;          % N*T
numThetaRowsTotal = Nvals * numThetaRowsLocal;     % N*M_i

% monomials in N*n_i variables of degree <= r_i
if degreeHi == 2
    numMonomialsTotal = (numStatesTotal + 2) .* (numStatesTotal + 1) / 2;
else
    numMonomialsTotal = arrayfun(@(nn) nchoosek(nn + degreeHi, degreeHi), ...
                                 numStatesTotal);
end

% same formulas evaluated at the aggregated quantities
numDecisionVarsMono = horizonTotal .* numStatesTotal .* numMonomialsTotal ...
                    + 1.5 .* numStatesTotal .* (numStatesTotal + 1) ...
                    + numStatesTotal.^2;

iterFactorMono = sqrt(2*numStatesTotal + horizonTotal) * log(1/epsi);

J_mono = iterFactorMono .* (numDecisionVarsMono.^3);

C_mono = numDecisionVarsMono .* (2*numStatesTotal + horizonTotal).^2 ...
       + numDecisionVarsMono.^2 ...
       + (horizonTotal .* numStatesTotal .* numMonomialsTotal ...
          + numStatesTotal .* (numStatesTotal + 1)/2) ...
         .* (numThetaRowsTotal .* numStatesTotal .* numMonomialsTotal);

%% ===================== Plots ============================================

% space complexity
plot_comparison(Nvals, C_comp, C_mono, ...
    '$\mathcal{C}_{\mathrm{comp}}$', '$\mathcal{C}_{\mathrm{mono}}$');

% time complexity
plot_comparison(Nvals, J_comp, J_mono, ...
    '$\mathcal{J}_{\mathrm{comp}}$', '$\mathcal{J}_{\mathrm{mono}}$');

%% ===================== Local function ===================================

function plot_comparison(Nvals, yComp, yMono, legComp, legMono)
    figure
    hold on; grid on; box on;
    loglog(Nvals, yComp, 'LineWidth', 2);
    valid = isfinite(yMono) & (yMono > 0);
    loglog(Nvals(valid), yMono(valid), 'LineWidth', 2);
    set(gca, 'TickLabelInterpreter', 'latex', 'XScale', 'log', 'YScale', 'log');
    xlabel('Number of Subsystems', 'Interpreter', 'latex');
    ylabel('Computational Cost', 'Interpreter', 'latex');
    legend({legComp, legMono}, 'Location', 'northwest', 'Interpreter', 'latex');
    xlim([1 2000]);
end