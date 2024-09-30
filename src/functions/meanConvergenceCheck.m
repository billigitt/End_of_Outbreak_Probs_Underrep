function [probFutureCases, risksBySimulation] = meanConvergenceCheck(GibbsSamples, numSamples)

% This function is simply used in the testing procedure to find the sampling values such that we
% get stable, non-changing risk estimates. The test is on a very simple epidemic, with only 2 cases.
% We can then copare the outputs to see if we have converged at the same estimates.

%% SI calculation

SI_mean = 15.3/7;
SI_sd = 9.3/7;

SI_scale = SI_sd^2/SI_mean;
SI_shape = SI_mean/SI_scale;


SI_discrete = zeros(50,1);
for k = 1:50
    intVals = [k-1:(1/10000):k+1];
    funcVals = (1 - abs(intVals - k)).*gampdf(intVals,SI_shape,SI_scale);
    SI_discrete(k) = trapz(intVals, funcVals);
end

SI_discrete(1) = SI_discrete(1) + (1-sum(SI_discrete));

%% Synthetic data generation

I_1 = 10;
rho = 0.5;
T = 5;
RpreERT = 0.5;
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

rhoAssumed = rho;% 0.6 0.7 0.8 0.9 0.95];
burnin = (GibbsSamples/10)-1;
thinning = 10;


for i = 1:length(rhoAssumed)

    outputGibbs = GibbsApproach([C; zeros(tEndCalc-T, 1)], SI_discrete, rhoAssumed(i), GibbsSamples, RpreERT, RERT, tERTArrival, tEndCalc, burnin, thinning);

    probFutureCases = outputGibbs.probMoreCases;

end

Tinf = 500;
experimentsPerSample = 1e3;

% simulationOutput = empiricalUnderrepEOOSmallMem(C, RpreERT, RERT, rhoAssumed(1), SI_discrete', numSamples, tERTArrival, tEndCalc, Tinf, experimentsPerSample);
simulationOutput = empiricalUnderrepEOOSmallMem(C, RpreERT, RERT, rhoAssumed(1), SI_discrete', numSamples, tERTArrival, tEndCalc, Tinf, experimentsPerSample);

risksBySimulation = simulationOutput.riskWithdrawERT;

%%
