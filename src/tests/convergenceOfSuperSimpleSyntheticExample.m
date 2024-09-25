% Script to find the sample values that yield stable probability estimates
% for the two different methods

%For this  specific example, we want to compare idx 3:end for both answers.

addpath('../functions')

GibbsSamples = 1e5;
simSamples = 1e1;

[probFutureCases1e5, risksBySimulation1e1] = meanConvergenceCheck(GibbsSamples, simSamples);

