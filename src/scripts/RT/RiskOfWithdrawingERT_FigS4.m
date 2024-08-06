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

number_missed = 1;
realisation_number = 1;

number_missed
realisation_number 

rng(realisation_number);


%% Consider the larger EVD dataset

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
ERTWithdrawalDay = daysact('5-apr-2018',  '24-jul-2018') + 1;


% Discretise a gamma distributed serial interval distribution using the method of
% Cori et al.;
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


%Now calculate the likelihood of R using the data up until the day before the arrival of the ERT assuming a constant R across this time period
RVals = [0.01:0.01:8];
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
likelihoodVals_Before_ERT = likelihoodVals./normfac; % Normalising the likelihood so that it is a valid pdf



%Now calculate the likelihood of R from ERT arrival until the day before
%ERT withdrawal
logL = zeros(length(RVals),1);

for rIt = 1:length(RVals)
    R = RVals(rIt);
    for t = ERTArrivalDay:(ERTWithdrawalDay - 1)
        k = incidence_data(t);
        lambda = 0;
        for s = 1:(t - 1)
            if s < (length(SI_discrete) + 1)
            lambda = lambda + R*incidence_data(t - s)*SI_discrete(s);
            end
        end
        logL(rIt) = logL(rIt) + k*log(lambda) - lambda - log(factorial(k));
    end
end

logL = logL - max(logL);
likelihoodVals = exp(logL);

normfac = trapz(RVals,likelihoodVals);
likelihoodVals_With_ERT = likelihoodVals./normfac; % Normalising the likelihood so that it is a valid pdf

likelihoodVals = likelihoodVals_Before_ERT;



% Now use the distributional estimates of R to construct the relative
% probability that, if a missed case occurred, that is was missed
% on each day of the outbreak

likelihoodOfMissedCase = zeros(1000,1);

for t = 2:(ERTArrivalDay - 1)
  
    for rIt = 1:length(RVals)
    R = RVals(rIt);
   for s = 1:(t - 1)
            if s < (length(SI_discrete) + 1)
                likelihoodOfMissedCase(t) = likelihoodOfMissedCase(t) + R*likelihoodVals_Before_ERT(rIt)*incidence_data(t - s)*SI_discrete(s);
            end
        end
    end
end

for t = ERTArrivalDay:(ERTWithdrawalDay - 1)
    
    for rIt = 1:length(RVals)
    R = RVals(rIt);
   for s = 1:(t - 1)
            if s < (length(SI_discrete) + 1)
                likelihoodOfMissedCase(t) = likelihoodOfMissedCase(t) + R*likelihoodVals_With_ERT(rIt)*incidence_data(t - s)*SI_discrete(s);
            end
        end
    end
end
likelihoodOfMissedCase = likelihoodOfMissedCase./sum(likelihoodOfMissedCase);

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


% Now calculate the risk of withdrawing the ERT each day from the day of
% the arrival of the ERT onwards - using the entire likelihood
% function for R (rather than a single value of R). Note that EOO (the end
% of outbreak probability) is 1 - Risk of withdrawing ERT.

% First assign missed cases

if number_missed > 0
    cumulativeProbs = cumsum(likelihoodOfMissedCase);
for caseNo = 1:number_missed


    r = rand();
    x = find(cumulativeProbs > r);
    dayMissed = x(1);

    incidence_data(dayMissed) = incidence_data(dayMissed) + 1;

end
end


% Calculate the likelihood of R using the data up until the day before
% the arrival of the ERT assuming a constant R across this time period, for
% the new "dataset" that includes the missed cases
RVals = [0.01:0.01:8];
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
likelihoodVals_Before_ERT = likelihoodVals./normfac; % Normalising the likelihood so that it is a valid pdf
likelihoodVals = likelihoodVals_Before_ERT;


% Now compute the risk of withdrawing the ERT 
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


%%% Now plot the figure (unlike in the manuscript, in this code this is only for the single
%%% (number_missed, realisation_number) pair specified at the beginning of the script)
ERTArrivalDay = daysact('5-apr-2018',  '8-may-2018') + 1;
ERTWithdrawalDay = daysact('5-apr-2018',  '24-jul-2018') + 1;

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