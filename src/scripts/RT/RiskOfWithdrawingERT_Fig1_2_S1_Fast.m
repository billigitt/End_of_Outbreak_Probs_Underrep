clear all
load('../../mats/varsFig1_2_S1.mat')

figure()
bar([1:61], [incidence_data(1:60); 0], 'BarWidth', 1)
hold on
plot([ERTArrivalDay ERTArrivalDay], [0 2], 'k--')
xticks([1:10:61])
% xticklabels(datesHere(1:10:61,:))
yticks([0 1 2])
box off
xlabel('Date')
ylabel('Number of cases')

figure()
bar([1:length(SI_discrete)], SI_discrete, 'BarWidth', 1)
xlim([0.5 50.5])
xticks([1 [5:5:50]])
box off

figure()
plot(RVals, likelihoodVals)
xlim([0 5])
box off

figure()
plot([0:56], 1 - combinedEOO(ERTArrivalDay:(ERTArrivalDay + 56)))
serial_first_case = datenum('15-may-2017');
serials = [serial_first_case:(serial_first_case+51)];
datesHere = datestr(serials);
xticks([0:10:50])
xticklabels(datesHere(1:10:51,:))
xlim([0 55])
box off
hold on
plot([48 48], [0 1], 'k--')