clc
clear
close all

numStatesLocal    = 3;
horizonLocal      = 17;
degreeHi          = 2;
numThetaRowsLocal = 9;
epsi              = 1e-8;

Nvals = 1:2000;

numMonomialsLocal = nchoosek(numStatesLocal + degreeHi, degreeHi);

numDecisionVarsLocal = horizonLocal*numStatesLocal*numMonomialsLocal ...
                     + 1.5*numStatesLocal*(numStatesLocal + 1) ...
                     + numStatesLocal^2;

iterFactorLocal = sqrt(2*numStatesLocal + horizonLocal) * log(1/epsi);

J_i_time = iterFactorLocal * (numDecisionVarsLocal^3);
J_comp   = Nvals .* J_i_time;

C_i_space = numDecisionVarsLocal*(2*numStatesLocal + horizonLocal)^2 ...
          + numDecisionVarsLocal^2 ...
          + (horizonLocal*numStatesLocal*numMonomialsLocal ...
             + numStatesLocal*(numStatesLocal + 1)/2) ...
            * (numThetaRowsLocal*numStatesLocal*numMonomialsLocal);
C_comp    = Nvals .* C_i_space;

numStatesTotal    = Nvals * numStatesLocal;
horizonTotal      = Nvals * horizonLocal;
numThetaRowsTotal = Nvals * numThetaRowsLocal;

if degreeHi == 2
    numMonomialsTotal = (numStatesTotal + 2) .* (numStatesTotal + 1) / 2;
else
    numMonomialsTotal = arrayfun(@(nn) nchoosek(nn + degreeHi, degreeHi), ...
                                 numStatesTotal);
end

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

plot_comparison(Nvals, C_comp, C_mono, ...
    '$\mathcal{C}_{\mathrm{comp}}$', '$\mathcal{C}_{\mathrm{mono}}$');

plot_comparison(Nvals, J_comp, J_mono, ...
    '$\mathcal{J}_{\mathrm{comp}}$', '$\mathcal{J}_{\mathrm{mono}}$');

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
