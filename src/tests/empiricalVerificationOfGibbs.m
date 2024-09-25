clear all
close all
clc

addpath('../functions')
rng(3)

% This script is the main body of the novel analysis that we conduct. We
% seek to use Gibbs sampling to account for under-reporting when
% calculating the end of outbreak probability.

%% SI calculation

SI_mean = 15.3/7;
SI_sd = 9.3/7;

SI_scale = SI_sd^2/SI_mean;
SI_shape = SI_mean/SI_scale;


SI_discrete = zeros(20,1); % I think life is made easier when calculating likelihoods 
% if the SI is atleast as long as tEndCalc
for k = 1:20
    intVals = [k-1:(1/10000):k+1];
    funcVals = (1 - abs(intVals - k)).*gampdf(intVals,SI_shape,SI_scale);
    SI_discrete(k) = trapz(intVals, funcVals);
end

%% Synthetic data generation

GibbsSamples = 1e5;
I_1 = 10;
rho = 0.5;
T = 5;
RpreERT = 1;
RERT = 0.5;
tERTArrival = floor(2*T/3);

I = zeros(T, 1);
I(1) = I_1;

for t = 2:(tERTArrival-1)

    I(t) = poissrnd(RpreERT*dot(I((t-1):-1:1), SI_discrete(1:(t-1))));

end

for t = tERTArrival:T

    I(t) = poissrnd(RERT*dot(I((t-1):-1:1), SI_discrete(1:(t-1))));

end

C = binornd(I, rho*ones(T, 1));
C(:) = 0;
C(1) = 2;
C(end) = 1;
C = [1 0 1]';
T = length(C);

tERTArrival = 3;

%

% RT original method (no under-reporting)

siCumulative = cumsum(SI_discrete);

[shape, rate] = posteriorRtEntireTimeSeries(C(1: (tERTArrival-1)), siCumulative);

scale = 1/rate;

riskTrueRCaseOnly = zeros(100, 1);
riskDistributionalRCaseOnly = zeros(100, 1);
gammaCaseOnly = zeros(100, 1);

for t = tERTArrival:1000

    if t <= (1+length(C))

        assumedCases = [C(1:(t-1)); zeros(100-t+1, 1)];

    else

        assumedCases = [C; zeros(100 - length(C), 1)];

    end

    gammaCaseOnly(t) = gamma_tMostBasic(assumedCases, SI_discrete, t);
    riskDistributionalRCaseOnly(t) = riskWithdrawERT(shape, rate, gammaCaseOnly(t));
    riskTrueRCaseOnly(t) = 1 - exp(-RpreERT*gammaCaseOnly(t));

end

tEndCalc = 15;

figure
yyaxis left
bar(C, 'BarWidth', 1)
hold on
ylabel('Reported cases')
yyaxis right
h(1) = plot(tERTArrival:tEndCalc, riskDistributionalRCaseOnly(tERTArrival:tEndCalc), 'g', 'LineWidth', 2, 'LineStyle', '-');
h(2) = plot(tERTArrival:tEndCalc, riskTrueRCaseOnly(tERTArrival:tEndCalc), 'b', 'LineWidth', 2, 'LineStyle', '-');

rhoAssumed = rho;% 0.6 0.7 0.8 0.9 0.95];
burnin = (GibbsSamples/10)-1;
thinning = 10;


% load('../mats/GibbsApproach1-Sensible-N1e4-CompareWithWrongEmpirical')



for i = 1:length(rhoAssumed)

    outputGibbs = GibbsApproach1Edit([C; zeros(tEndCalc-T, 1)], SI_discrete, rhoAssumed(i), GibbsSamples, RpreERT, RERT, tERTArrival, tEndCalc, burnin, thinning);

    probFutureCases = outputGibbs.probMoreCases;

    h(i+2) = plot(tERTArrival:tEndCalc, probFutureCases(tERTArrival:tEndCalc), 'color', i*[1 1 1]/(1+length(rhoAssumed)), 'LineWidth', 2, 'LineStyle', '-');

end

%%

numSamples = 2e5;
Tinf = 200;
experimentsPerSample = 1e3;

% simulationOutput = empiricalUnderrepEOOSmallMem(C, RpreERT, RERT, rhoAssumed(1), SI_discrete', numSamples, tERTArrival, tEndCalc, Tinf, experimentsPerSample);
tic
simulationOutput = empiricalUnderrepEOOSmallMem(C, RpreERT, RERT, rhoAssumed(1), SI_discrete', numSamples, tERTArrival, tEndCalc, Tinf, experimentsPerSample);
toc
risksBySimulation = simulationOutput.riskWithdrawERT;

%
h(4) = plot(risksBySimulation);

legend(h([1 2 3 4]), 'RT et. al. (distributional R)', 'RT et. al. (true R)', 'Gibbs rho=0.5', 'risks by sim')%, ...
%     'Gibbs rho=0.6', 'Gibbs rho=0.7', 'Gibbs rho=0.8', 'Gibbs rho=0.9',...
%     'Gibbs rho=0.95', 'Gibbs rho=0.99')
xlim([0 tEndCalc])
ylabel("Probability of future"+newline+"cases after day t")
xlabel('Day t')
xline(tERTArrival, '--')
