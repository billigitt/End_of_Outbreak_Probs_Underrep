% Script to find the sample values that yield stable probability estimates
% for the two different methods

%For this  specific example, we want to compare idx 3:end for both answers.

addpath('../functions')

GibbsSamples = 1e3;
simSamples = 2e5;

[~, risksBySimulation2e5] = meanConvergenceCheck(GibbsSamples, simSamples);

simSamples = 2e4;
[~, risksBySimulation2e4] = meanConvergenceCheck(GibbsSamples, simSamples);

