clc
clear all
close all

%% This code generates the results shown in Figure S3 of the manuscript:
%% Thompson et al. Using real-time modelling to inform the 2017 Ebola outbreak response in DR Congo. Nature Communications, 2024.

%% All code was written in MATLAB, compatible with version R2022a.

%% ©2024 Robin Thompson <robin.thompson@maths.ox.ac.uk>

incidence_data = zeros(1000,1);
incidence_data(1) = 1;

incidence_data(daysact('5-apr-2018',  '8-apr-2018') + 1) = incidence_data(daysact('5-apr-2018',  '8-apr-2018') + 1) + 1;
incidence_data(daysact('5-apr-2018',  '12-apr-2018') + 1) = incidence_data(daysact('5-apr-2018',  '12-apr-2018') + 1) + 3;
incidence_data(daysact('5-apr-2018',  '13-apr-2018') + 1) = incidence_data(daysact('5-apr-2018',  '13-apr-2018') + 1) + 1;
incidence_data(daysact('5-apr-2018',  '18-apr-2018') + 1) = incidence_data(daysact('5-apr-2018',  '18-apr-2018') + 1) + 1;
incidence_data(daysact('5-apr-2018',  '19-apr-2018') + 1) = incidence_data(daysact('5-apr-2018',  '19-apr-2018') + 1) + 1;
incidence_data(daysact('5-apr-2018',  '20-apr-2018') + 1) = incidence_data(daysact('5-apr-2018',  '20-apr-2018') + 1) + 1;
incidence_data(daysact('5-apr-2018',  '21-apr-2018') + 1) = incidence_data(daysact('5-apr-2018',  '21-apr-2018') + 1) + 1;
incidence_data(daysact('5-apr-2018',  '23-apr-2018') + 1) = incidence_data(daysact('5-apr-2018',  '23-apr-2018') + 1) + 1;
incidence_data(daysact('5-apr-2018',  '24-apr-2018') + 1) = incidence_data(daysact('5-apr-2018',  '24-apr-2018') + 1) + 1;
incidence_data(daysact('5-apr-2018',  '25-apr-2018') + 1) = incidence_data(daysact('5-apr-2018',  '25-apr-2018') + 1) + 1;
incidence_data(daysact('5-apr-2018',  '27-apr-2018') + 1) = incidence_data(daysact('5-apr-2018',  '27-apr-2018') + 1) + 1;
incidence_data(daysact('5-apr-2018',  '1-may-2018') + 1) = incidence_data(daysact('5-apr-2018',  '1-may-2018') + 1) + 2;
incidence_data(daysact('5-apr-2018',  '3-may-2018') + 1) = incidence_data(daysact('5-apr-2018',  '3-may-2018') + 1) + 2;
incidence_data(daysact('5-apr-2018',  '4-may-2018') + 1) = incidence_data(daysact('5-apr-2018',  '4-may-2018') + 1) + 6;
incidence_data(daysact('5-apr-2018',  '5-may-2018') + 1) = incidence_data(daysact('5-apr-2018',  '5-may-2018') + 1) + 2;
incidence_data(daysact('5-apr-2018',  '6-may-2018') + 1) = incidence_data(daysact('5-apr-2018',  '6-may-2018') + 1) + 1;
incidence_data(daysact('5-apr-2018',  '7-may-2018') + 1) = incidence_data(daysact('5-apr-2018',  '7-may-2018') + 1) + 1;
incidence_data(daysact('5-apr-2018',  '8-may-2018') + 1) = incidence_data(daysact('5-apr-2018',  '8-may-2018') + 1) + 1;
incidence_data(daysact('5-apr-2018',  '10-may-2018') + 1) = incidence_data(daysact('5-apr-2018',  '10-may-2018') + 1) + 3;
incidence_data(daysact('5-apr-2018',  '12-may-2018') + 1) = incidence_data(daysact('5-apr-2018',  '12-may-2018') + 1) + 1;
incidence_data(daysact('5-apr-2018',  '13-may-2018') + 1) = incidence_data(daysact('5-apr-2018',  '13-may-2018') + 1) + 2;
incidence_data(daysact('5-apr-2018',  '14-may-2018') + 1) = incidence_data(daysact('5-apr-2018',  '14-may-2018') + 1) + 2;
incidence_data(daysact('5-apr-2018',  '15-may-2018') + 1) = incidence_data(daysact('5-apr-2018',  '15-may-2018') + 1) + 3;
incidence_data(daysact('5-apr-2018',  '16-may-2018') + 1) = incidence_data(daysact('5-apr-2018',  '16-may-2018') + 1) + 3;
incidence_data(daysact('5-apr-2018',  '18-may-2018') + 1) = incidence_data(daysact('5-apr-2018',  '18-may-2018') + 1) + 2;
incidence_data(daysact('5-apr-2018',  '19-may-2018') + 1) = incidence_data(daysact('5-apr-2018',  '19-may-2018') + 1) + 2;
incidence_data(daysact('5-apr-2018',  '20-may-2018') + 1) = incidence_data(daysact('5-apr-2018',  '20-may-2018') + 1) + 3;
incidence_data(daysact('5-apr-2018',  '21-may-2018') + 1) = incidence_data(daysact('5-apr-2018',  '21-may-2018') + 1) + 2;
incidence_data(daysact('5-apr-2018',  '28-may-2018') + 1) = incidence_data(daysact('5-apr-2018',  '28-may-2018') + 1) + 1;
incidence_data(daysact('5-apr-2018',  '2-jun-2018') + 1) = incidence_data(daysact('5-apr-2018',  '2-jun-2018') + 1) + 1;

ERTArrivalDay = daysact('5-apr-2018',  '8-may-2018') + 1;


% Results are plotted both for the outbreak-specific serial interval (SI =
% 1) and the original serial interval (SI = 2)

for SI = 1

if SI == 1
SI_mean = 19.46;
SI_sd = 6.08;
end
if SI == 2
SI_mean = 15.3;
SI_sd = 9.3;
end

% Discretise a gamma distributed serial interval distribution using the method of
% Cori et al.;
SI_scale = SI_sd^2/SI_mean;
SI_shape = SI_mean/SI_scale;


SI_discrete = zeros(100,1);
for k = 1:100
    intVals = [k-1:(1/10000):k+1];
    funcVals = (1 - abs(intVals - k)).*gampdf(intVals,SI_shape,SI_scale);
    SI_discrete(k) = trapz(intVals, funcVals);
end
SI_discrete = SI_discrete./sum(SI_discrete);

%Now calculate the likelihood of R using the data up until the day before the arrival of the ERT assuming a constant R across this time period
RVals = [0.01:0.01:10];
logL = zeros(length(RVals),1);

for rIt = 1:length(RVals)
    R = RVals(rIt);
    for t = 2:(ERTArrivalDay - 1)
        k = incidence_data(t);
        lambda = 0;
        for s = 1:(t - 1)
            lambda = lambda + R*incidence_data(t - s)*SI_discrete(s);
        end
        logL(rIt) = logL(rIt) + k*log(lambda) - lambda - log(factorial(k));
    end
end

logL = logL - max(logL);
likelihoodVals = exp(logL);

normfac = trapz(RVals,likelihoodVals);
likelihoodVals = likelihoodVals./normfac; % Normalising the likelihood so that it is a valid pdf


% Now calculate the risk of withdrawing the ERT each day from the day of
% the arrival of the ERT onwards - using the entire likelihood
% function for R (rather than a single value of R). Note that EOO (the end
% of outbreak probability) is 1 - Risk of withdrawing ERT.
logEOO = zeros(ERTArrivalDay + 150,length(RVals));
EOO = zeros(ERTArrivalDay + 150,length(RVals));
for RIt = 1:length(RVals)
    R = RVals(RIt)
    for t = ERTArrivalDay:(ERTArrivalDay + 150)
        
        incidence_data_here = incidence_data;
        incidence_data_here((t+1):length(incidence_data_here)) = 0;

        logEOO(t,RIt) = 0;
        for j = t:(t + length(SI_discrete))
            for s = 1:(j-1)
                if s < length(SI_discrete + 0.5)
                    logEOO(t,RIt) = logEOO(t,RIt) - R*incidence_data_here(j-s)*SI_discrete(s);
                end
            end
        end

        EOO(t,RIt) = exp(logEOO(t, RIt));

    end
end
combinedEOO = zeros(ERTArrivalDay + 150,1);
for t = ERTArrivalDay:(ERTArrivalDay + 150)
    funcV = zeros(length(RVals),1);
    for i = 1:length(RVals)
        funcV(i) = EOO(t,i)*likelihoodVals(i);
    end
    combinedEOO(t) = trapz(RVals, funcV);
end


%%% Now plot the figures
serial_first_case = datenum('5-apr-2018');
serials = [serial_first_case:(serial_first_case+71)];
datesHere = datestr(serials);

figure()
bar([1:71], [incidence_data(1:70); 0], 'BarWidth', 1)
hold on
plot([ERTArrivalDay ERTArrivalDay], [0 6], 'k--')
xticks([1:10:71])
xticklabels(datesHere(1:10:71,:))
yticks([0:6])
box off
xlabel('Date')
ylabel('Number of cases')
xlim([0 61])

figure()
bar([1:length(SI_discrete)], SI_discrete, 'BarWidth', 1)
xlim([0.5 50.5])
xticks([1 [5:5:50]])
box off

figure()
plot(RVals, likelihoodVals)
xlim([0 8])
box off

figure(4)
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

end