clc
clear all
close all

%% This code relates to the results in Figure S4 of the manuscript:
%% Thompson et al. Using real-time modelling to inform the 2017 Ebola outbreak response in DR Congo. Nature Communications, 2024.

%% Specifically, to generate the results in Fig 4B, this code should be run for number_missed = 0,1,2,...,5
%% For each value of number_missed, the code should be run for realisation_number = 1,2,...,1000
%% The results shown in Fig 4B are the average values of 1 - combinedEOO across all 1000 realisations for each fixed value of number_missed

%% All code was written in MATLAB, compatible with version R2022a.

%% ©2024 Robin Thompson <robin.thompson@maths.ox.ac.uk>

clc
clear all
close all

load('../../mats/varsS4.mat')

figure()
bar([1:1000], likelihoodOfMissedCase, 'BarWidth', 1)
hold on
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

figure()
hold on
plot([0:100], 1 - combinedEOO(ERTArrivalDay:(ERTArrivalDay + 100)))
serial_first_case = datenum('8-may-2018');
serials = [serial_first_case:(serial_first_case+101)];
datesHere = datestr(serials);
xticks([0:20:100])
xticklabels(datesHere(1:20:101,:))
xlim([0 100])
yticks([0:0.2:1])
box off
hold on
last_case = daysact('8-may-2018',  '24-jul-2018');
plot([last_case last_case], [0 1], 'k--')