clear all
close all
clc

addpath('../functions')
rng(2)

% This script is the main body of the novel analysis that we conduct. We
% seek to use Gibbs sampling to account for under-reporting when
% calculating the end of outbreak probability.

%% SI calculation (copied from RT's github)

SI_mean = 15.3;
SI_sd = 9.3;

SI_scale = SI_sd^2/SI_mean;
SI_shape = SI_mean/SI_scale;


SI_discrete = zeros(100,1);
for k = 1:100
    intVals = [k-1:(1/10000):k+1];
    funcVals = (1 - abs(intVals - k)).*gampdf(intVals,SI_shape,SI_scale);
    SI_discrete(k) = trapz(intVals, funcVals);
end

%% Synthetic data generation

GibbsSamples = 1e2;
I_1 = 10;
rho = 0.5;
T = 50;
RpreERT = 3;
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

%

% RT original method (no under-reporting)

siCumulative = cumsum(SI_discrete);

[shape, rate] = posteriorRtEntireTimeSeries(C(1: (tERTArrival-1)), siCumulative);

scale = 1/rate;

riskTrueRCaseOnly = zeros(1000, 1);
riskDistributionalRCaseOnly = zeros(1000, 1);
gammaCaseOnly = zeros(1000, 1);

for t = tERTArrival:1000
    
    if t <= length(C)

        assumedCases = [C(1:t); zeros(1000-t, 1)];

    else

        assumedCases = [C; zeros(1000 - length(C), 1)];

    end

    gammaCaseOnly(t) = gamma_tMostBasic(assumedCases, SI_discrete, t);
    riskDistributionalRCaseOnly(t) = riskWithdrawERT(shape, rate, gammaCaseOnly(t));
    riskTrueRCaseOnly(t) = 1 - exp(-RpreERT*gammaCaseOnly(t));

end


tEndCalc = 1e2;

figure
yyaxis left
bar(C, 'BarWidth', 1)
hold on
ylabel('Reported cases')
yyaxis right
h(1) = plot(tERTArrival:tEndCalc, riskDistributionalRCaseOnly(tERTArrival:tEndCalc), 'g', 'LineWidth', 2, 'LineStyle', '-');
h(2) = plot(tERTArrival:tEndCalc, riskTrueRCaseOnly(tERTArrival:tEndCalc), 'b', 'LineWidth', 2, 'LineStyle', '-');

rhoAssumed = rho;% 0.6 0.7 0.8 0.9 0.95];
burnin = round(GibbsSamples/10);


for i = 1:length(rhoAssumed)

tic 
profile on
outputGibbs = GibbsApproach1([C; zeros(tEndCalc-T, 1)], SI_discrete, rhoAssumed(i), GibbsSamples, RpreERT, RERT, tERTArrival, tEndCalc, burnin);
toc

probFutureCases = outputGibbs.probMoreCases;

h(i+2) = plot(tERTArrival:tEndCalc, probFutureCases(tERTArrival:tEndCalc), 'color', i*[1 1 1]/(1+length(rhoAssumed)), 'LineWidth', 2, 'LineStyle', '-');

end

numSamples = 1e3;
Tinf = 300;
experimentsPerSample = 1e3;

tic
risksBySimulation = empiricalUnderrepEOO(C, RpreERT, RERT, rhoAssumed(1), SI_discrete', numSamples, tERTArrival, tEndCalc, Tinf, experimentsPerSample);
toc

h(4) = plot(risksBySimulation);

%

legend(h([1 2 3 4]), 'RT et. al. (distributional R)', 'RT et. al. (true R)', 'Gibbs rho=0.5', 'risks by sim')%, ...
%     'Gibbs rho=0.6', 'Gibbs rho=0.7', 'Gibbs rho=0.8', 'Gibbs rho=0.9',...
%     'Gibbs rho=0.95', 'Gibbs rho=0.99')
xlim([0 tEndCalc])
ylabel("Probability of future"+newline+"cases after day t")
xlabel('Day t')
xline(tERTArrival, '--')