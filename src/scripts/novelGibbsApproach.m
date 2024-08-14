clear all
close all
clc

addpath('../functions')
rng(1)

% This script is the main body of the novel analysis that we conduct. We
% seek to use Gibbs sampling to account for under-reporting when
% calculating the end of outbreak probability.

%% SI calculation (copied from RT's github)

SI_mean = 7.3;
SI_sd = 6.3;

SI_scale = SI_sd^2/SI_mean;
SI_shape = SI_mean/SI_scale;


SI_discrete = zeros(100,1);
for k = 1:100
    intVals = [k-1:(1/10000):k+1];
    funcVals = (1 - abs(intVals - k)).*gampdf(intVals,SI_shape,SI_scale);
    SI_discrete(k) = trapz(intVals, funcVals);
end

%% Synthetic data generation

GibbsSamples = 1e1;
I_1 = 5;
rho = 0.3;
T = 50;
RpreERT = 2;
RERT = 0.5;
tERTArrival = floor(T/3);

I = zeros(T, 1);
I(1) = I_1;

for t = 2:(tERTArrival-1)

    I(t) = poissrnd(RpreERT*dot(I((t-1):-1:1), SI_discrete(1:(t-1))));

end

for t = tERTArrival:T

    I(t) = poissrnd(RERT*dot(I((t-1):-1:1), SI_discrete(1:(t-1))));

end

C = binornd(I, rho*ones(T, 1));

%%

tEndCalc = 1e2;
tic
profile on
probFutureCases = GibbsApproach1([C; zeros(tEndCalc-T, 1)], SI_discrete, rho, GibbsSamples, RpreERT, RERT, tERTArrival, tEndCalc);
toc

%% RT original method (no under-reporting)

siCumulative = cumsum(SI_discrete);

[shape, rate] = posteriorRtEntireTimeSeries(C(1: (tERTArrival-1)), siCumulative);

scale = 1/rate;

riskCaseOnly = zeros(1000, 1);
gammaCaseOnly = zeros(1000, 1);

for t = tERTArrival:1000
    
    if t <= length(C)

        assumedCases = [C(1:t); zeros(1000-t, 1)];

    else

        assumedCases = [C; zeros(1000 - length(C), 1)];

    end

    gammaCaseOnly(t) = gamma_tMostBasic(assumedCases, SI_discrete, t);
    riskCaseOnly(t) = riskWithdrawERT(shape, rate, gammaCaseOnly(t));
    
end

%%

figure
bar(I, 'BarWidth', 1)
hold on
plot(probFutureCases, 'r')
plot(riskCaseOnly, 'b')