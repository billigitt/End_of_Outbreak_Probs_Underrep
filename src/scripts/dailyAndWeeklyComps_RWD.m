% Cleaning
clear all
clc
close all

set(0,'DefaultFigureWindowStyle', 'normal')
set(0,'DefaultTextInterpreter', 'none')
set(0,'defaultAxesTickLabelInterpreter','none')
set(0, 'defaultLegendInterpreter','none')
set(0, 'defaultaxesfontsize', 22)
set(0, 'defaultlinelinewidth', 1.5)

set(0, 'DefaultAxesFontName', 'aakar');
set(0, 'DefaultTextFontName', 'aakar');
set(0, 'defaultUicontrolFontName', 'aakar');
set(0, 'defaultUitableFontName', 'aakar');
set(0, 'defaultUipanelFontName', 'aakar');

set(groot, 'defaultAxesTickLabelInterpreter','tex');
set(groot, 'defaultLegendInterpreter','tex');
set(0, 'DefaultTextInterpreter', 'tex')

colourMat = [0.9 0.6 0;... %orange 1
    0.35 0.7 0.9;... %sky blue 2
    0 0.6 0.5;... %blueish green 3
    0.9290 0.6940 0.1250;... %yellow 4
    0 0.44 0.7;... %blue 5
    0.8 0.4 0;... %red 6
    0.4940 0.1840 0.5560]; % purple 7

%Risk(t) is in this file are defined to be the probability of no more cases
% after (and including) week t, given data (full weeks) from week the first week to week
% t-1.

addpath('../functions')

rng(1)

% Daily serial interval calc

SI_mean_daily = 15.3;
SI_sd_daily = 9.3;

SI_scale_daily = SI_sd_daily^2/SI_mean_daily;
SI_shape_daily = SI_mean_daily/SI_scale_daily;


wDaily = zeros(100,1);
for k = 1:100
    intVals = [k-1:(1/10000):k+1];
    funcVals = (1 - abs(intVals - k)).*gampdf(intVals,SI_shape_daily,SI_scale_daily);
    wDaily(k) = trapz(intVals, funcVals);
end

wDaily(1) = wDaily(1) + (1-sum(wDaily));

% Weekly serial interval calc

P = 7;
SI_mean_weekly = 15.3/P;
SI_sd_weekly = 9.3/P;

SI_scale_weekly = SI_sd_weekly^2/SI_mean_weekly;
SI_shape_weekly = SI_mean_weekly/SI_scale_weekly;

wWeekly = zeros(100,1);
for k = 1:100
    intVals = linspace((k-1), (k+1), 1000);
    funcVals = (1 - abs(intVals - (k))).*gampdf(intVals,SI_shape_weekly,SI_scale_weekly);
    wWeekly(k) = trapz(intVals, funcVals);
end

wWeekly(1) = wWeekly(1) + (1-sum(wWeekly));

% Daily incidence import, with Rt = 2 (pre ERT), Rt = 0.5 (during ERT)

%The first day of the epidemiological week is 2nd April 2018. And so we
%want this date to correspond to index 1. Since datenum('02-apr-2018') =
%737152, when "datenumming" each date, we then minus 737151 so that the
%indices are correctly alligned.

totalTime = 2.1e2;
incidenceData = zeros(totalTime, 1);
idx1s = datenum([2018 4 5; 2018 4 8; 2018 4 13; 2018 4 18;...
    2018 4 19; 2018 4 20; 2018 4 21; 2018 4 23; 2018 4 24; 2018 4 25;...
    2018 4 27; 2018 5 6; 2018 5 7; 2018 5 8; 2018 5 12; 2018 5 28; ...
    2018 6 2]) - 737151;
idx2s = datenum([2018 5 1; 2018 5 3; 2018 5 5; 2018 5 13; 2018 5 14; ...
    2018 5 18; 2018 5 19; 2018 5 21]) - 737151;
idx3s = datenum([2018 4 12; 2018 5 10; 2018 5 15; 2018 5 16; 2018 5 20])...
    - 737151;
idx6s = datenum([2018 5 4]) - 737151;
incidenceData(idx1s') = 1;
incidenceData(idx2s') = 2;
incidenceData(idx3s') = 3;
incidenceData(idx6s') = 6;

idxERTDeployed = datenum([2018 5 8]) - 737151;
idxERTWithdrawn = datenum([2018 7 24]) - 737151;
incidenceBeforeERT = incidenceData(1:(idxERTDeployed-1));

idxRiskPlotEnd = datenum([2018 9 15]) - 737151;

% Rename as cases

rho = 0.5;
CDaily = incidenceData;

% Weekly case calc (summing from daily cases). Also generate schematic.
% FIG 1

CDailyMat = reshape(CDaily, 7, totalTime/7);
CWeekly = sum(CDailyMat)';

weeks = [datetime(2018, 4, 5, 12, 0, 0): caldays(7): datetime(2018, 10, 25, 12, 0, 0)] ;%This will take the first day of each epi week up to length of C
figure
bar(weeks(1:18), CWeekly(1:18), 'BarWidth',1, 'LineStyle','none')
hold on
xtickformat('dd/MM')
xticks(datetime(2018, 4, 2) : caldays(7) : datetime(2018, 10, 22))
xtickangle(45)

fill([datetime(2018, 5, 8) datetime(2018, 7, 24) datetime(2018, 7, 24)...
    datetime(2018,5, 8)], [0 0 30 30], 'r', 'FaceAlpha',0.5, 'LineStyle','none')
text(datetime(2018, 6, 15), 22, "Actual"+newline+"ERT period", 'HorizontalAlignment','center', 'Color','r', 'FontSize',25, 'FontName','aakar')
% text(datetime(2018, 7, 24), 22, "ERT"+newline+"withdrawal date", 'HorizontalAlignment','center', 'Color','r', 'FontSize', 20)
text(datetime(2018, 6, 30), 12, "Withdraw ERT?", 'Color', 'k', 'HorizontalAlignment','center', 'FontSize',25, 'FontName','aakar')
ylim([0 20])

for i = 1:13

    x = [0.39925 + 0.514*(i-1)/13, 0.39925 + 0.514*(i-1)/13];
    y = [0.5 0.31];

    annotation('arrow', x, y)

end

ylabel('Reported incidence')
xlabel('Date (dd/mm)')



ax = gca;
labels = string(ax.XAxis.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
ax.XAxis.TickLabels = labels; % set
title(" "+newline+" ")
box off
set(gcf,'Position',[100 100 1150 600])
set(gcf, 'color', 'none')
%% Calculate Rt posteriors (pre ERT and during ERT) (and using weekly and daily data)

%% Daily calculation - currently we do this such that we technically have more informationfor the daily inference than the weekly.
CDailyPreERT = CDaily(4:(7*floor((idxERTDeployed-1)/7))); %4 because 4th index correponds to first case
CDailyERT = CDaily(7*ceil(idxERTDeployed/7)+1:(idxERTWithdrawn-1));

siCumulative = cumsum(wDaily);

[shapeBeforeERTDaily, rateBeforeERTDaily] = posteriorRtEntireTimeSeries(CDailyPreERT, siCumulative);
scaleBeforeERTDaily = 1/rateBeforeERTDaily;
meanRBeforeERTDaily = shapeBeforeERTDaily*scaleBeforeERTDaily;

shapeDuringERT = 1 + sum(CDailyERT);
rateDuringERT = gamma_tMostBasic(incidenceData(1:(idxERTWithdrawn-1)), wDaily, idxERTDeployed);
scaleDuringERT = 1/rateDuringERT;
meanRDuringERTDaily = shapeDuringERT*scaleDuringERT;

%% Weekly calculation
%idxERTDeployed = 37, and so since 37/7 = 5.29, weeks 1 to 5 are pre-ERT,
%and weeks 7 onwards are during ERT. Week 6is a hybrid week so cannot be
%deemed either. idxERTWithdrawn/7 = 16.29 and so we only infer R from weeks
%7 to 16 for during ERT.
CWeeklyPreERT = CWeekly(1:5);
CWeeklyERT = CWeekly(7:16);

siCumulativeWeekly = cumsum(wWeekly);

[shapeBeforeERTWeekly, rateBeforeERTWeekly] = posteriorRtEntireTimeSeries(CWeeklyPreERT, siCumulativeWeekly);
scaleBeforeERTWeekly = 1/rateBeforeERTWeekly;
meanRBeforeERTWeekly = shapeBeforeERTWeekly*scaleBeforeERTWeekly;

shapeDuringERT = 1 + sum(CWeeklyERT);
rateDuringERT = gamma_tMostBasic(CWeekly, wWeekly, ceil(idxERTDeployed/7)); %check that the third arg is correct here
scaleDuringERT = 1/rateDuringERT;
meanRDuringERTWeekly = shapeDuringERT*scaleDuringERT;

%% Calculate EOO with daily data (using end pts + daily SI)

TinfDaily = 7*30;
tEndAnalysis = 7*30 + 1; %The +1 means that you estimate the probability of 0s from the next day after
% a full week, which will be comparable to weekly data

gammaDailyCalc = zeros(tEndAnalysis, 1);
riskDailyCalc = gammaDailyCalc;


for t = (7*ceil(idxERTDeployed/7)+1):tEndAnalysis %Starting at day 43.

    if t <= (1+length(CDaily))

        assumedCases = [CDaily(1:(t-1)); zeros(TinfDaily-t+1, 1)]; %first non-zero part (for first loop entry) should be of length 42 since first 7 weeks form the data

    else

        assumedCases = [CDaily; zeros(TinfDaily - length(CDaily), 1)];

    end

    gammaDailyCalc(t) = gamma_tMostBasic(assumedCases, wDaily, t);
    riskDailyCalc(t) = 1 - exp(-meanRBeforeERTDaily*gammaDailyCalc(t)); %Calculates the probability of atleast one more case for day t onwards, given days 1:t-1

end

 % we want all contributions to be mod 1
% as these will be the ones that use the data from the  full week

%% Calculate EOO with weekly data (using weekly SI)

TinfWeekly = 50;
tEndAnalysisWeekly = 31;
tERTArrivalWeekly = ceil(idxERTDeployed/7)+1;

gammaWeeklyCalc = zeros(tEndAnalysisWeekly, 1);
riskWeeklyCalc = gammaWeeklyCalc;

casesCheckWeek = [];

for t = tERTArrivalWeekly:tEndAnalysisWeekly

    if t <= (1+length(CWeekly))

        assumedCases = [CWeekly(1:t-1); zeros(TinfWeekly-t+1, 1)]; %This should be of length 6 for the non-zero part on first loop.

    else

        assumedCases = [CWeekly; zeros(TinfWeekly - length(CWeekly), 1)];

    end

    gammaWeeklyCalc(t) = gamma_tMostBasic(assumedCases, wWeekly, t);
    riskWeeklyCalc(t) = 1 - exp(-meanRBeforeERTWeekly*gammaWeeklyCalc(t));

end

%% Plot to compare


days = [datetime(2018, 4, 2): caldays(1): datetime(2018, 10, 28)];
weeks = [datetime(2018, 4, 2): caldays(7): datetime(2018, 10, 29)];
figure
bar(days(1:65), CDaily(1:65), 'LineStyle', 'none', 'BarWidth', 1)
xtickformat('dd/MM')
xticks(datetime(2018, 4, 2) : caldays(7) : datetime(2018, 10, 22))
xtickangle(45)
ax = gca;
labels = string(ax.XAxis.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
ax.XAxis.TickLabels = labels; % set
xlabel('Date (dd/mm)')
ylabel('Reported incidence')
set(gcf,'Position',[100 100 1150 600])
set(gcf, 'color', 'none')
box off


figure
% subplot(2, 1, 1)
% bar(days, CDaily)
% xtickformat('dd-MMM')
% xticks(datetime(2018, 4, 2) : caldays(7*3) : datetime(2018, 10, 22))
% xlabel('Days')
% ylabel('Daily C_t')
% subplot(2, 2, 2)
% bar(weeks, CWeekly)
% xlabel('Weeks')
% ylabel('Weekly C_t')
% xtickformat('dd-MMM')
% xticks(datetime(2018, 4, 2) : caldays(7*3) : datetime(2018, 10, 22))
subplot(1, 2, 1)
RR = linspace(0, 8, 1e3);
h(1) = plot(RR, gampdf(RR, shapeBeforeERTWeekly, scaleBeforeERTWeekly), 'color', colourMat(5, :));
hold on
xline(meanRBeforeERTWeekly, 'color', colourMat(5, :), 'LineWidth', 2)
h(2) = plot(RR, gampdf(RR, shapeBeforeERTDaily, scaleBeforeERTDaily), 'color', colourMat(6, :));
xline(meanRBeforeERTDaily, 'color', colourMat(6, :), 'LineWidth', 2)
xlabel('\itR\rm_{pre-ERT}')
ylabel('Probability density')
legend(h([2 1]), 'Daily data', 'Weekly data')
xlim([0 8])
box off
subplot(1, 2, 2)
%riskDailyCalc is indexed by day, and we want day 43 to end in jumps of 7
h(1) = plot(weeks(ceil(idxERTDeployed/7)+1:31), riskDailyCalc(43:7:end), 'color', colourMat(6, :)); %Take off first value because the first point for daily is with 0 data, where as weekly has a whole week.
hold on
%riskWeeklyCalc is indexed by week, and we want week 7 to end.
h(2) = plot(weeks(ceil(idxERTDeployed/7)+1:31), riskWeeklyCalc(7:end), 'color', colourMat(5, :), 'LineStyle', '--'); 
legend(h([1 2]), 'Daily data', 'Weekly data')
xlabel('Date (dd/mm)')
ylabel("Probability of"+newline+"future cases")
xtickformat('dd/MM')
xticks(weeks(ceil(idxERTDeployed/7)+1:31))
xlim([weeks(ceil(idxERTDeployed/7)+1) weeks(31)])
xtickangle(45)
ax = gca;
labels = string(ax.XAxis.TickLabels); % extract
logicals = (rem(1:length(labels), 3)~=1);
labels(logicals) = nan; % remove every other one
labels(logicals) = nan;
ax.XAxis.TickLabels = labels; % set
box off
set(gcf,'Position',[100 100 1150 600])
set(gcf, 'color', 'none')

%% Using weekly data, account for under-reporting.

rho = [(0.3:0.1:0.8) 1];
GibbsSamples = 1e5;
burnin = GibbsSamples/10;
thinning = 10;
probMoreCasesWithRho = zeros(tEndAnalysisWeekly, length(rho));
statGeweke = zeros(tEndAnalysisWeekly, length(rho));

%Un-comment these lines when computing risks again

% for i = 1:length(rho)
% 
%     disp(i)
%     output = GibbsApproach1Edit(CWeekly, wWeekly, rho(i), GibbsSamples, meanRBeforeERTWeekly, meanRDuringERTWeekly, ...
%         tERTArrivalWeekly, tEndAnalysisWeekly, burnin, thinning);
% 
%     probMoreCasesWithRho(:, i) = output(1).probMoreCases;
%     statGeweke(:, i) = output(1).pGeweke;
% end

load('../mats/fig3Gibbs1e5')

weekSafe = zeros(7, 1);

figure
subplot(1, 2, 1)
hold on
for i = 1:6

    if i == 4

        h(i) = plot(weeks(ceil(idxERTDeployed/7)+1:30), probMoreCasesWithRho(ceil(idxERTDeployed/7)+1:end, i), 'color', ...
            ((8-i)/8)*colourMat(2, :), 'Marker', 'x');

    else

        h(i) = plot(weeks(ceil(idxERTDeployed/7)+1:30), probMoreCasesWithRho(ceil(idxERTDeployed/7)+1:end, i), 'color', ...
            ((8-i)/8)*colourMat(2, :));
    end
    weekSafe(i) = find(probMoreCasesWithRho(6:end, i)<=0.01, 1);

end

h(7) = plot(weeks(ceil(idxERTDeployed/7)+1:30), probMoreCasesWithRho(ceil(idxERTDeployed/7)+1:end, 7), 'k--');
weekSafe(7) = find(probMoreCasesWithRho(6:end, 7)<=0.01, 1);

xtickformat('dd-MMM')
xticks(datetime(2018, 4, 2) : caldays(7*3) : datetime(2018, 10, 22))
ylabel("Probability of"+newline+"future cases")
xlabel('Date')
h(8) = xline(datetime(2018, 7, 24), 'r', 'LineWidth', 2);
legend(h(1:8), '\fontsize{15}\rho = 0.3', '\fontsize{15}\rho = 0.4', ...
    '\fontsize{15}\rho = 0.5', '\fontsize{15}\rho = 0.6', '\fontsize{15}\rho = 0.7',...
    '\fontsize{15}\rho = 0.8', '\fontsize{15}\rho = 1', ...
    "\fontsize{15}24/07 (Actual"+newline+"ERT with-"+newline+"drawal date)")
xlim([weeks(ceil(idxERTDeployed/7)+1) weeks(30)])
xtickformat('dd/MM')
subplot(1, 2, 2)
bar(rho, weekSafe, 'FaceColor',colourMat(2, :))
xlabel('Reporting probability, \rho')
ylabel('Theoretical ERT period')
set(gcf,'Position',[100 100 1350 600])
set(gcf, 'color', 'none')
box off
xlim([0.25 1.05])


%% FIG S1

load('../mats/convergenceOfGibbs.mat')
load('../mats/convergenceOfSims.mat')
%%
figure
yyaxis left
bar([[1 0 1] zeros(1, length(probFutureCases1e5)-3)], 'LineStyle', 'none', 'BarWidth', 1)
ylabel('Reported cases')
yticks([0 1])
hold on
yyaxis right
h(1) = plot([[nan; nan]; probFutureCases1e5(3:end)], 'color', colourMat(2, :), 'LineWidth', 2);
h(2) = plot([[nan; nan]; risksBySimulation2e5(3:end-2)], 'k--', 'LineWidth', 2);
xlabel('Weeks')
ylabel("Probability of"+newline+"future cases")
ylim([0 1])
xtickangle(0)
xlim([0 13])

ax = gca;
ax.YAxis(2).Color = 'k';
box off

set(gcf,'Position',[100 100 700 600])
set(gcf, 'color', 'none')

legend(h([1 2]), 'Novel method', 'Simulation method')

disp("max error for Gibbs between 1e6 and 1e5 is "+max(abs(probFutureCases1e5 - probFutureCases1e6)))
disp("max error for Gibbs between 1e5 and 1e4 is "+max(abs(probFutureCases1e5 - probFutureCases1e4)))

disp("max error for Simulations between 2e4 and 2e5 is "+max(abs(risksBySimulation2e5 - risksBySimulation2e4)))

disp("max error in plots is"+max(abs(probFutureCases1e5(3:end) - risksBySimulation2e5(3:end-2))))

% This suggests that we need to choose Gibbs = 1e5 and Sims = 2e4 (does 2e4
% look wierd in the manuscript?) I am now thinking that an absolute error
% of less than 0.005 is fine (low relative errors will require really high
% number of 
% samples)
%NB: used rho = 0.5, R = 0.5, need to redo Gibbs because we need to
%estimate risk for slightly longer.