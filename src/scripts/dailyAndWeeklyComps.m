%% Cleaning
clear all
clc
close all


addpath('../functions')

rng(1)

%% Daily serial interval calc

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

%% Weekly serial interval calc

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

%% Daily incidence calc, with Rt = 2 (pre ERT), Rt = 0.5 (during ERT)

I_1 = 10;
T = 7*11;
I = [I_1; zeros(T-1, 1)];
tERTArrivalDaily = 7*4+1;
tERTArrivalWeekly = 4;
RpreERT = 2;
RERT = 0.5;

for t = 2:T

    R = RpreERT + (RERT - RpreERT)*(t>=tERTArrivalDaily);
    I(t) = renewalEqn(I(1:(t-1)), wDaily', R);

end

%% Convert to cases

rho = 0.5;
CDaily = binornd(I, rho*ones(T, 1));

%% Weekly case calc (summing from daily cases)

CDailyMat = reshape(CDaily, 7, T/7);
CWeekly = sum(CDailyMat)';

%% Calculate EOO with daily data (using end pts + daily SI)

TinfDaily = 7*30;
tEndAnalysis = 7*20 + 1;

gammaDailyCalc = zeros(tEndAnalysis, 1);
riskDailyCalc = gammaDailyCalc;

casesCheckDay = [];

for t = tERTArrivalDaily:tEndAnalysis

    if t <= (1+length(CDaily))

        assumedCases = [CDaily(1:(t-1)); zeros(TinfDaily-t+1, 1)];

    else

        assumedCases = [CDaily; zeros(TinfDaily - length(CDaily), 1)];

    end

    if (rem(t, 7) == 1)

        casesCheckDay = [casesCheckDay sum(assumedCases)];

    end

    gammaDailyCalc(t) = gamma_tMostBasic(assumedCases, wDaily, t);
    riskDailyCalc(t) = 1 - exp(-RpreERT*gammaDailyCalc(t));

end

riskDailyCalc = riskDailyCalc(1:7:tEndAnalysis); % we want all contributions to be mod 1
% as these will be the ones that use the data from the entire full week

%% Calculate EOO with weekly data (using weekly SI)

TinfWeekly = 30;
tEndAnalysisWeekly = 20;

gammaWeeklyCalc = zeros(tEndAnalysisWeekly, 1);
riskWeeklyCalc = gammaWeeklyCalc;

casesCheckWeek = [];

for t = tERTArrivalWeekly:tEndAnalysisWeekly

    if t <= length(CWeekly)

        assumedCases = [CWeekly(1:t); zeros(TinfWeekly-t+1, 1)];

    else

        assumedCases = [CWeekly; zeros(TinfWeekly - length(CWeekly), 1)];

    end

    casesCheckWeek = [casesCheckWeek sum(assumedCases)];

    gammaWeeklyCalc(t) = gamma_tMostBasic(assumedCases, wWeekly, t+1);
    riskWeeklyCalc(t) = 1 - exp(-RpreERT*gammaWeeklyCalc(t));

end
%% Plot to compare

figure
subplot(2, 2, 1)
bar(CDaily)
xlabel('Days')
ylabel('Daily case count')
subplot(2, 2, 2)
bar(CWeekly)
xlabel('Weeks')
ylabel('Weekly case count')
subplot(2, 2, 3)
xlabel('R')
ylabel('pdf (for weekly and daily data)')
subplot(2, 2, 4)
h(1) = plot(riskDailyCalc(2:end), 'b'); %Take off first value because the first point for daily is with 0 data, where as weekly has a whole week.
hold on
h(2) = plot(riskWeeklyCalc, 'r', 'LineStyle', '--'); 
legend(h([1 2]), 'Daily Risk Calc', 'Weekly Risk Calc')
xlabel('After seeing t weeks of data')
ylabel('Probability of seeing no more cases after (non inclusive) week t')

figure
plot(casesCheckWeek, 'b')
hold on
plot(casesCheckDay, 'r', 'LineStyle','--')
ylabel('Probability of no more cases from week t onwards')
xlabel('Week t')

%% Using weekly data, account for under-reporting.