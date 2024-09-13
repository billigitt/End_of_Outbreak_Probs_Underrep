close all
clc
clear all

load('../mats/GibbsApproach1-SuperSimple-N1e5-CompareWithWrongEmpirical')



tt = 2;

for t = 1:4

    if t+1 >= tt

        figure
        h(1) = histogram(simulationOutput(t+1).cellOfMatrices(tt, :), 'BinWidth', 1, 'Normalization','probability');
        hold on
        h(2) = histogram(outputGibbs(t).incidenceStore(tt, :), 'BinWidth', 1, 'Normalization','probability');

        legend(h([1 2]), 'sims', 'gibbs')

    end

end
