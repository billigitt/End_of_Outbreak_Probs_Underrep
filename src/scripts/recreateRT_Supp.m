clear all
close all

%This file is similar to the recreateRT.m file, except that we now look at
%recreating the supplementary figures from the same manuscript. In
%particular, we wish to recreate the results from Figures S4a and S4b.

%incidence in 2018: 27-mar(1), 18-apr(1), 24-apr(2), 26-apr(1), 11-may(1), 30-apr(1), 02-may(1)
%ERT arrival in 2018: 15-may

%day 1 is 05-apr-2018, and datenum('05-apr-2018') = 737155. So we can index
%all dates to be datenum('dd-mm-2018')-737154. We will make the entire
%incidence 1000 in length (very long).

addpath('../functions')

load('../mats/varsS4.mat', 'likelihoodOfMissedCase', 'ERTArrivalDay',...
    'ERTWithdrawalDay', 'likelihoodVals_With_ERT') %load RT's calculations from S4

load('../mats/varsS3.mat', 'likelihoodVals', 'RVals') %load RT's calculations from S3


tic

totalTime = 3.1e2;
incidenceData = zeros(totalTime, 1);
idx1s = datenum([2018 4 5; 2018 4 8; 2018 4 13; 2018 4 18;...
    2018 4 19; 2018 4 20; 2018 4 21; 2018 4 23; 2018 4 24; 2018 4 25;...
    2018 4 27; 2018 5 6; 2018 5 7; 2018 5 8; 2018 5 12; 2018 5 28; ...
    2018 6 2]) - 737154;
idx2s = datenum([2018 5 1; 2018 5 3; 2018 5 5; 2018 5 13; 2018 5 14; ...
    2018 5 18; 2018 5 19; 2018 5 21]) - 737154;
idx3s = datenum([2018 4 12; 2018 5 10; 2018 5 15; 2018 5 16; 2018 5 20])...
    - 737154;
idx6s = datenum([2018 5 4]) - 737154;
incidenceData(idx1s') = 1;
incidenceData(idx2s') = 2;
incidenceData(idx3s') = 3;
incidenceData(idx6s') = 6;

idxERTDeployed = datenum([2018 5 8]) - 737154;
idxERTWithdrawn = datenum([2018 7 24]) - 737154;
incidenceBeforeERT = incidenceData(1:(idxERTDeployed-1));

idxRiskPlotEnd = datenum([2018 9 15]) - 737154;

%% SI calculation (copied from RT's github)

%NB: The serial interval that we use here is from the maximum likelihood
%calculation (see figure S3C)
SI_mean = 19.46;
SI_sd = 6.08;

SI_scale = SI_sd^2/SI_mean;
SI_shape = SI_mean/SI_scale;


SI_discrete = zeros(100,1);
for k = 1:100
    intVals = [k-1:(1/10000):k+1];
    funcVals = (1 - abs(intVals - k)).*gampdf(intVals,SI_shape,SI_scale);
    SI_discrete(k) = trapz(intVals, funcVals);
end
SI_discrete = SI_discrete./sum(SI_discrete);


%%

siCumulative = cumsum(SI_discrete);

[shapeBeforeERT, rateBeforeERT] = posteriorRtEntireTimeSeries(incidenceBeforeERT, siCumulative);
scaleBeforeERT = 1/rateBeforeERT;
meanRBeforeERT = shapeBeforeERT*scaleBeforeERT;

shapeDuringERT = 1 + sum(incidenceData(idxERTDeployed:(idxERTWithdrawn-1)));
rateDuringERT = gamma_tMostBasic(incidenceData(1:(idxERTWithdrawn-1)), SI_discrete, idxERTDeployed);
scaleDuringERT = 1/rateDuringERT;
meanRDuringERT = shapeDuringERT*scaleDuringERT;

numberOfMs = 5;
numberOfSamples = 1e3;

risks = zeros(totalTime, numberOfSamples, numberOfMs);
riskAggregates = zeros(totalTime, numberOfMs);
gamma = zeros(totalTime, numberOfSamples, numberOfMs);

relativeFOI = zeros(idxERTWithdrawn-1, 1);
for t = 1:(idxERTDeployed-1)

    relativeFOI(t) = forceOfInfection(incidenceData, SI_discrete, meanRBeforeERT, t);

end

for t = idxERTDeployed:(idxERTWithdrawn-1)

    relativeFOI(t) = forceOfInfection(incidenceData, SI_discrete, meanRDuringERT, t);

end

%Now turn into probability distribution by normalizing
relativeFOI = relativeFOI/sum(relativeFOI);

samplesOfHiddenInfections = zeros(numberOfSamples, totalTime, numberOfMs);

for m = 1:numberOfMs

    samplesOfHiddenInfections(:, 1:(idxERTWithdrawn-1), m) = mnrnd(m, relativeFOI, numberOfSamples);

end

samplesOfHiddenInfections = permute(samplesOfHiddenInfections, [2 1 3]);

samplesOfTrueInfections = samplesOfHiddenInfections + repmat(incidenceData, [1 numberOfSamples numberOfMs]);

%Calculate shape and scale for each sampled true incidence

[shapeBeforeERTSample, rateBeforeERTSample] = posteriorRtEntireTimeSeries3D(samplesOfTrueInfections(1:(idxERTDeployed-1), :, :), siCumulative);

for t = idxERTDeployed:totalTime

    gamma(t, :, :) = gamma_tMostBasic3D(samplesOfTrueInfections, SI_discrete, t);
    risks(t, :, :) = riskWithdrawERT(shapeBeforeERTSample, rateBeforeERTSample, squeeze(gamma(t, :, :)));
    riskAggregates(t, :) = squeeze(mean(risks(t, :, :), 2));

end

%% Plots

figure
bar(incidenceData, 'BarWidth', 1, 'LineStyle','-')
xlim([0 60])

figure
x = linspace(0,10,1000);
plot(x, gampdf(x, shapeBeforeERT, scaleBeforeERT), 'LineWidth', 2)
hold on
plot(RVals, likelihoodVals, '--', 'LineWidth', 2)


RVals = [0.01:0.01:8];

figure
x = linspace(0,10,1000);
plot(x, gampdf(x, shapeDuringERT, scaleDuringERT), 'LineWidth', 2)
hold on
plot(RVals, likelihoodVals_With_ERT, '--', 'LineWidth', 2)

figure
bar([1:1000], likelihoodOfMissedCase, 'BarWidth', 1, 'FaceAlpha',0.5)
hold on
bar((1:length(relativeFOI))', relativeFOI, 'BarWidth', 1, 'FaceAlpha', 0.5, 'LineStyle','none')
plot([ERTArrivalDay ERTArrivalDay], [0 6], 'k--')
plot([ERTWithdrawalDay ERTWithdrawalDay], [0 6], 'k--')
xlim([0 300])
ylim([0 0.05])

serial_first_case = datenum('5-apr-2018');
serials = [serial_first_case:(serial_first_case+121)];
datesHere = datestr(serials);
xticks([0:20:120])
xticklabels(datesHere(1:20:121,:))
xlim([0 120])
box off
hold on

figure
hold on
for m = 1:numberOfMs

    h(m) = plot(1:totalTime, riskAggregates(:, m));
    
end

legend(h(1:numberOfMs), 'm=1', 'm=2', 'm=3', 'm=4', 'm=5')

% 
% x = linspace(0, 5, 1000);
% xx = (1:1000);
% 
% figure
% plot(x, gampdf(x, shapeBeforeERT, scaleBeforeERT), 'LineWidth', 2)
% ylabel('Prob density')
% xlabel('Reproduction number (R)')
% 
% figure
% h(1) = plot(737154+xx(idxERTDeployed: idxRiskPlotEnd), risk(xx(idxERTDeployed: idxRiskPlotEnd)), 'LineWidth', 2);
% hold on
% h(2) = plot(737154+xx(idxERTDeployed: idxRiskPlotEnd), 1.-combinedEOO(xx(idxERTDeployed: idxRiskPlotEnd)), 'LineWidth', 2, 'LineStyle', '--');
% datetick('x', 19, "keepticks")
% ylabel('Risk of withdrawing ERT')
% xlabel('Date (DD/MM)')
% legend(h([1 2]), 'ZOG', 'RT')
% 
% figure
% h(1) = plot(xx, gamma);
% hold on
% h(2) = plot(xx, gamma_tRT, '--');
% legend(h([1 2]), 'ZOG', 'RT')
% ylabel('Gamma')
% xlabel('t')
% 
toc
