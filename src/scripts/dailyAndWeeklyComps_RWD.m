% ------------
% We consider a disease incidence time series dataset from an EVD outbreak that took place in Équateur Province, DRC, in 2018.
% In total, 54 cases (40 confirmed cases and 14 probable cases) occurred between 5th April and 2nd June 2018. 
% The Emergency Response Team (ERT) was deployed on 8th May 2018 and was withdrawn on 24th July 2018. 
% The overall goal of our analysis is to determine the risk of withdrawing the ERT (i.e., the probability of future cases) at the beginning of each week following the deployment of the ERT.
% ------------

% Cleaning
clear all
clc
% close all

%% Set figure specifications
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



%% Make available required functions and package files
addpath('../functions')
addpath('../packages/mcmcstat-master') % use mcmcstat for geweke calculation

%% Seed the random number generator
rng(1)

%% Daily serial interval calc

% Specify mean and standard deviation of gamma-distributed serial interval (SI)
SI_mean_daily = 15.3;
SI_sd_daily = 9.3;

% Parameterise gamma distributed SI in terms of scale and shape parameters
SI_scale_daily = SI_sd_daily^2/SI_mean_daily;
SI_shape_daily = SI_mean_daily/SI_scale_daily;

% Discretise the SI (daily timesteps)
wDaily = zeros(100,1);
for k = 1:100
    intVals = [k-1:(1/10000):k+1];
    funcVals = (1 - abs(intVals - k)).*gampdf(intVals,SI_shape_daily,SI_scale_daily);
    wDaily(k) = trapz(intVals, funcVals);
end

% Numerical correction to first entry to ensure wDaily sums to 1
wDaily(1) = wDaily(1) + (1-sum(wDaily));

%% Weekly serial interval calc

% Specify mean and standard deviation of gamma-distributed serial interval (SI)
P = 7;% Conversion factor from days to weeks of 7 days
SI_mean_weekly = 15.3/P;
SI_sd_weekly = 9.3/P;

% Parameterise gamme distributed SI in terms of scale and shape parameters
SI_scale_weekly = SI_sd_weekly^2/SI_mean_weekly;
SI_shape_weekly = SI_mean_weekly/SI_scale_weekly;

% Discretise the SI (weekly timesteps)
wWeekly = zeros(100,1);
for k = 1:100
    intVals = linspace((k-1), (k+1), 1000);
    funcVals = (1 - abs(intVals - (k))).*gampdf(intVals,SI_shape_weekly,SI_scale_weekly);
    wWeekly(k) = trapz(intVals, funcVals);
end

% Numerical correction to first entry to ensure wDaily sums to 1
wWeekly(1) = wWeekly(1) + (1-sum(wWeekly));

%% Data set up

% ------------
% Disease incidence data: 
%  - EVD outbreak that took place in Équateur Province, DRC, in 2018.
%  - 54 cases (40 confirmed cases and 14 probable cases) occurred between 5th April and 2nd June 2018. 
%  - ERT deployed on 8th May 2018 and was withdrawn on 24th July 2018. 
% ------------

% ------------
% The first day of the epidemiological week was 2nd April 2018. And so we
% want this date to correspond to index 1. Since datenum('02-apr-2018') =
% 737152, when "datenumming" each date, we then minus 737151 so that the
% indices are correctly alligned.
% ------------

% Set up time horizon for analysis & initialise storage vector for
% incidence counts
totalTime = 2.1e2;
incidenceData = zeros(totalTime, 1);

% Set up array entry access values for case counts (1 case, 2 cases, 3 cases, 6 cases) that occurred on specified dates
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

% Set up array entry access values corresponding to ERT periods
idxERTDeployed = datenum([2018 5 8]) - 737151;
idxERTWithdrawn = datenum([2018 7 24]) - 737151;
ERTPeriod = idxERTWithdrawn - idxERTDeployed; % Specify length of ERT period

% Assign incidence for time period before ERT was deployed to vector
incidenceBeforeERT = incidenceData(1:(idxERTDeployed-1));

% Case reporting proportion parameter
rho = 0.5;

%% Generate FIG 1 - visulisation of the ERT withdrawl time schematic

%rename incidenceData as CDaily
CDaily = incidenceData;

% Weekly case calc (summing from daily cases). Also generate schematic.
%% Generate FIG 1 - visulisation of the ERT withdrawl time schematic

% Weekly case calc (summing from daily cases).
CDailyMat = reshape(CDaily, 7, totalTime/7);
CWeekly = sum(CDailyMat)';

% Plot generation
weeks = [datetime(2018, 4, 5, 12, 0, 0): caldays(7): datetime(2018, 10, 25, 12, 0, 0)] ;%This will take the first day of each epi week up to length of C
figure
fill([datetime(2018, 5, 8) datetime(2018, 7, 24) datetime(2018, 7, 24)...
    datetime(2018,5, 8)], [0 0 30 30], 'r', 'FaceAlpha',0.25, 'LineStyle','none')
hold on
xtickformat('dd/MM')
xticks(datetime(2018, 4, 2) : caldays(7) : datetime(2018, 10, 22))
xtickangle(45)
bar(weeks(1:18), CWeekly(1:18), 'BarWidth',1, 'LineStyle','none')

text(datetime(2018, 6, 15), 22, "Actual ERT period", 'HorizontalAlignment','center', 'Color','r', 'FontSize',25, 'FontName','aakar')
text(datetime(2018, 6, 30), 12, "Withdraw ERT?", 'Color', 'k', 'HorizontalAlignment','center', 'FontSize',25, 'FontName','aakar')
ylim([0 20])

for i = 1:13

    x = [0.39925 + 0.514*(i-1)/13, 0.39925 + 0.514*(i-1)/13];
    y = [0.5 0.31];

    annotation('arrow', x, y)

end

% Specify axis labels 
ylabel('Weekly reported cases')
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

% Daily calculation 
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

% Weekly R calculation
%idxERTDeployed = 37, and so since 37/7 = 5.29, weeks 1 to 5 are pre-ERT,
%and weeks 7 onwards are during ERT. Week 6 is a hybrid week so cannot be
%deemed either. idxERTWithdrawn/7 = 16.29 and so we only infer R from weeks
%7 to 16 for during ERT.
CWeeklyPreERT = CWeekly(1:5);
CWeeklyERT = CWeekly(7:16);

siCumulativeWeekly = cumsum(wWeekly);

[shapeBeforeERTWeekly, rateBeforeERTWeekly] = posteriorRtEntireTimeSeries(CWeeklyPreERT, siCumulativeWeekly);
scaleBeforeERTWeekly = 1/rateBeforeERTWeekly;
meanRBeforeERTWeekly = shapeBeforeERTWeekly*scaleBeforeERTWeekly;

shapeDuringERT = 1 + sum(CWeeklyERT);
rateDuringERT = gamma_tMostBasic(CWeekly, wWeekly, ceil(idxERTDeployed/7)+1); %3rd argument is 7, and corresponds to the week from which we multiply the poisson likelihoods together
scaleDuringERT = 1/rateDuringERT;
meanRDuringERTWeekly = shapeDuringERT*scaleDuringERT;

%% Calculate EOO with daily data (using end pts + daily SI)

TinfDaily = 7*30;
tEndAnalysis = 7*30 + 1; 
    % The +1 means that you estimate the probability of 0s from the next day after
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

%% Generate FIG S1 - Daily numbers of reported cases in the 2018 EVD outbreak in Équateur Province, DRC.
days = [datetime(2018, 4, 2): caldays(1): datetime(2018, 10, 28)];
weeks = [datetime(2018, 4, 2): caldays(7): datetime(2018, 10, 29)];

figure
bar(days(1:65), CDaily(1:65), 'LineStyle', 'none', 'BarWidth', 1)
xtickformat('dd/MM')
xticks(datetime(2018, 4, 2) : caldays(7) : datetime(2018, 10, 22))
xtickangle(45)
xlabel('Date (dd/mm)')
ylabel('Daily reported cases')
set(gcf,'Position',[100 100 1150 600])
set(gcf, 'color', 'none')
box off

%% Generate FIG 2 - Probability of future cases estimated from either weekly or daily disease incidence time series data, assuming perfect case reporting.
figure
subplot(1, 2, 1)
RR = linspace(0, 8, 1e3);
h(1) = plot(RR, gampdf(RR, shapeBeforeERTWeekly, scaleBeforeERTWeekly), ...
    'color', 'k', 'LineStyle', '--', 'LineWidth', 2);
hold on
xline(meanRBeforeERTWeekly, 'color', 'k', 'LineWidth', 2, 'LineStyle', '--')
h(2) = plot(RR, gampdf(RR, shapeBeforeERTDaily, scaleBeforeERTDaily), 'color', colourMat(5, :), 'LineWidth', 2);
xline(meanRBeforeERTDaily, 'color', colourMat(5, :), 'LineWidth', 2)
xlabel('Reproduction number (\itR\rm)')
ylabel('Probability density')
legend(h([2 1]), 'Daily data', 'Weekly data')
xlim([0 8])
box off
subplot(1, 2, 2)
%riskDailyCalc is indexed by day, and we want day 43 to end in jumps of 7
h(1) = plot(weeks(ceil(idxERTDeployed/7)+1:31), riskDailyCalc(43:7:end), 'color', colourMat(5, :), 'LineWidth', 2); %Take off first value because the first point for daily is with 0 data, where as weekly has a whole week.
hold on
%riskWeeklyCalc is indexed by week, and we want week 7 to end.
h(2) = plot(weeks(ceil(idxERTDeployed/7)+1:31), riskWeeklyCalc(7:end), 'color', 'k', 'LineStyle', '--', 'LineWidth', 2); 
xline(datetime(2018, 7, 24), 'r', 'LineWidth', 2)
legend(h([1 2]), 'Daily data', 'Weekly data')
xlabel('Date (dd/mm)')
ylabel("Probability of future cases")
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

rho = 0.3:0.1:1;
GibbsSamples = 1e5; 
burnin = GibbsSamples/10;
thinning = 10; 
probMoreCasesWithRho = zeros(tEndAnalysisWeekly, length(rho));
statGewekeRho = zeros(tEndAnalysisWeekly, length(rho));

probGibbsSample = zeros(tEndAnalysisWeekly, GibbsSamples, length(rho));

%UN-COMMENT FOLLOWING LINES WHEN COMPUTING P(T) CURVE AGAIN

for i = 1:length(rho)

    disp(i)
    output = GibbsApproach(CWeekly, wWeekly, rho(i), GibbsSamples, meanRBeforeERTWeekly, meanRDuringERTWeekly, ...
        tERTArrivalWeekly, tEndAnalysisWeekly, burnin, thinning);

    probMoreCasesWithRho(:, i) = output(1).probMoreCases;
    probGibbsSample(:, :, i) = output(1).probGibbsSample;
    statGewekeRho(:, i) = output(1).pGeweke;
end

%%

%save('../mats/fig3Gibbs1e5RhoNew.mat') %set YOURFILENAME to fig3Gibbs1e5RhoNew
% or similar, and remove future loads.

%load('../mats/fig3Gibbs1e5RhoNew')

% Burn-in with 10,000 did not initially work so we use burn-in = 20,000 and find Geweke tests are passed in all cases

zGeweke = zeros(31, 8);
pGewekeBuiltIn = zeros(31, 8);

burnin = 2e4;

for j = 1:length(rho)

    for jj = 7:31

        [zGeweke(jj, j), pGewekeBuiltIn(jj, j)] = geweke(probGibbsSample(jj, burnin + 1:thinning:end, j)');

    end

    probMoreCasesWithRho(:, j) = mean(probGibbsSample(:, burnin + 1:thinning:end, j), 2);

end

%When choosing this burn-in, all the Geweke stats exceed 0.05, meaning we
%can accept the null hypothesis that the first 10% of the chain (sans burn
%in) has the same mean as the last 50% of the chain (sans burn-in). Thsi
%suggests that the chain has converged.

%% Generate FIG 3 - Estimated probability of future cases accounting for case under-reporting.

weekSafe = zeros(8, 1);

figure
subplot(1, 2, 1)
hold on
for i = 1:7

    if i == 4

        h(i) = plot(weeks(ceil(idxERTDeployed/7)+1:30), probMoreCasesWithRho(ceil(idxERTDeployed/7)+1:end-1, i), 'color', ...
            ((9-i)/9)*colourMat(2, :), 'Marker', 'x');

    else

        h(i) = plot(weeks(ceil(idxERTDeployed/7)+1:30), probMoreCasesWithRho(ceil(idxERTDeployed/7)+1:end-1, i), 'color', ...
            ((9-i)/9)*colourMat(2, :));
    end
    weekSafe(i) = find(probMoreCasesWithRho(7:end, i)<=0.05, 1)-1/7;

end

h(8) = plot(weeks(ceil(idxERTDeployed/7)+1:30), probMoreCasesWithRho(ceil(idxERTDeployed/7)+1:end-1, end), 'k--');
weekSafe(8) = find(probMoreCasesWithRho(7:end, 8)<=0.05, 1)-1/7;

xtickformat('dd/MM')
xticks(datetime(2018, 4, 2) : caldays(7*3) : datetime(2018, 10, 22))
ylabel("Probability of future cases")
xlabel('Date (dd/mm)')
h(9) = xline(datetime(2018, 7, 24), 'r', 'LineWidth', 2);
legend(h(1:9), '\fontsize{15}\rho = 0.3', '\fontsize{15}\rho = 0.4', ...
    '\fontsize{15}\rho = 0.5', '\fontsize{15}\rho = 0.6', '\fontsize{15}\rho = 0.7',...
    '\fontsize{15}\rho = 0.8', '\fontsize{15}\rho = 0.9', '\fontsize{15}\rho = 1', ...
    "\fontsize{15}Actual ERT"+newline+"withdrawal date")
xlim([weeks(ceil(idxERTDeployed/7)+1) weeks(30)])
xtickformat('dd/MM')
subplot(1, 2, 2)
bar(rho, weekSafe, 'FaceColor',colourMat(2, :))
hold on
yline(ERTPeriod/7, 'r', 'LineWidth', 2)
xlabel('Reporting probability (\rho)')
ylabel('Theoretical ERT period (weeks)')
set(gcf,'Position',[100 100 1550 600])
set(gcf, 'color', 'none')
box off
xlim([0.25 1.05])
xticks(0.3:0.1:1)
xticklabels({'0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1'})


%% Generate FIG S2 - Comparison of results from the Gibbs sampling method (as used in the main text) and an alternative approach involving repeated model simulation. 
load('../mats/convergenceOfGibbs.mat')
load('../mats/convergenceOfSims.mat')
%
figure
yyaxis left
bar([[1 0 1] zeros(1, length(probFutureCases1e5)-3)], 'LineStyle', 'none', 'BarWidth', 1)
ylabel('Reported cases')
yticks([0 1])
hold on
yyaxis right
h(1) = plot(0.5:12.5, [[nan; nan; nan]; probFutureCases1e5(4:end)], 'color', colourMat(2, :), 'LineWidth', 2);
h(2) = plot(0.5:12.5, [[nan; nan; nan]; risksBySimulation2e5(4:end-2)], 'k--', 'LineWidth', 2);
xlabel('Time (weeks)')
ylabel("Probability of future cases")
ylim([0 1])
xtickangle(0)
xlim([0 13])

ax = gca;
ax.YAxis(2).Color = 'k';
box off

set(gcf,'Position',[100 100 700 600])
set(gcf, 'color', 'none')

legend(h([1 2]), 'Gibbs sampling', 'Simulation method')

disp("max error for Gibbs between 1e6 and 1e5 is "+max(abs(probFutureCases1e5 - probFutureCases1e6)))
disp("max error for Gibbs between 1e5 and 1e4 is "+max(abs(probFutureCases1e5 - probFutureCases1e4)))

disp("max error for Simulations between 2e4 and 2e5 is "+max(abs(risksBySimulation2e5 - risksBySimulation2e4)))

disp("max error in plots is"+max(abs(probFutureCases1e5(3:end) - risksBySimulation2e5(3:end-2))))

% This suggests that we need to choose Gibbs = 1e5 and Sims = 2e4.

%% Generate FIG 4 - Accounting for uncertainty in the case reporting probability, 𝝆, when estimating the probability of future cases.

rhoPMF = [0.1 0.15 0.25 0.25 0.15 0.1];
weightedProbFutureCases = rhoPMF*probMoreCasesWithRho(1:end-1, 3:end)';

figure
subplot(1, 2, 1)
bar(0.5:0.1:1, rhoPMF, 'BarWidth', 1)
ylim([0 0.275])
xlim([0.425 1.075])
xlabel('Reporting probability (\rho)')
ylabel('Probability')
xticks(0.5:0.1:1)
box off

subplot(1, 2, 2)
h(1) = plot(weeks(ceil(idxERTDeployed/7)+1:30), weightedProbFutureCases(7:end));
hold on
h(2) = plot(weeks(ceil(idxERTDeployed/7)+1:30), probMoreCasesWithRho(ceil(idxERTDeployed/7)+1:end-1, end), 'k--');
xlabel('Date (dd/mm) ')
ylabel('Probability of future cases')
xtickangle(45)
xtickformat('dd/MM')
xticks(datetime(2018, 4, 2) : caldays(7*3) : datetime(2018, 10, 22))

box off

set(gcf,'Position',[100 100 1150 600])
set(gcf, 'color', 'none')

xlim([weeks(ceil(idxERTDeployed/7)+1) weeks(30)])
xline(datetime(2018, 7, 24), 'r', 'LineWidth', 2)
legend(h([1 2]), 'Uncertain \rho', '\rho = 1')

