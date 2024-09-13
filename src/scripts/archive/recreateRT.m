clear all
close all
%incidence in 2017: 27-mar(1), 18-apr(1), 24-apr(2), 26-apr(1), 11-may(1), 30-apr(1), 02-may(1)
%ERT arrival in 2017: 15-may

%day 1 is 27-mar-2017, and datenum('27-mar-2017') = 736781. So we can index
%all dates to be datenum('dd-mm-2017')-736780. We will make the entire
%incidence 1000 in length (very long).

addpath('../../functions')

load('../../mats/varsFig1_2_S1_ZOGEdit.mat', 'combinedEOO', 'gamma_tRT') %load RT's risk calculation

incidenceData = zeros(1000, 1);
idx1s = datenum([2017 3 27; 2017 4 18; 2017 4 26; 2017 4 30; 2017 5 2; 2017 5 11]) - 736780; 
idx2s = datenum([2017 4 24]) - 736780;

incidenceData(idx1s') = 1;
incidenceData(idx2s') = 2;

idxERTDeployed = datenum([2017 5 15]) - 736780;
incidenceBeforeERT = incidenceData(1:(idxERTDeployed-1));

idxRiskPlotEnd = datenum([2017 7 15]) - 736780;

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


%%

siCumulative = cumsum(SI_discrete);

[shape, rate] = posteriorRtEntireTimeSeries(incidenceBeforeERT, siCumulative);

scale = 1/rate;

risk = zeros(1000, 1);
gamma = zeros(1000, 1);

for t = idxERTDeployed:1000
    
    gamma(t) = gamma_tMostBasic(incidenceData, SI_discrete, t);
    risk(t) = riskWithdrawERT(shape, rate, gamma(t));
    
end


%% Plots

x = linspace(0, 5, 1000);
xx = (1:1000);

figure
plot(x, gampdf(x, shape, scale), 'LineWidth', 2)
ylabel('Prob density')
xlabel('Reproduction number (R)')

figure
h(1) = plot(736780+xx(idxERTDeployed: idxRiskPlotEnd), risk(xx(idxERTDeployed: idxRiskPlotEnd)), 'LineWidth', 2);
hold on
h(2) = plot(736780+xx(idxERTDeployed: idxRiskPlotEnd), 1.-combinedEOO(xx(idxERTDeployed: idxRiskPlotEnd)), 'LineWidth', 2, 'LineStyle', '--');
datetick('x', 19, "keepticks")
ylabel('Risk of withdrawing ERT')
xlabel('Date (DD/MM)')
legend(h([1 2]), 'ZOG', 'RT')

figure
h(1) = plot(xx, gamma);
hold on
h(2) = plot(xx, gamma_tRT, '--');
legend(h([1 2]), 'ZOG', 'RT')
ylabel('Gamma')
xlabel('t')

% As we cab see, the plots match exactly between the two scripts,
% confirming that our code is good!

