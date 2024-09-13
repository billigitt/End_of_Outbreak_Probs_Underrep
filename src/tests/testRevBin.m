close all

numTrials = 1e6;

trials = randi([0 100], numTrials, 1);

successes = binornd(trials, 0.3*ones(numTrials, 1));

idxMatch = (successes == 1);

figure
histogram(trials(idxMatch), 'Normalization','probability', 'BinWidth',1)
hold on
histogram(revBin(1, 0.3, 1e4), 'Normalization','probability', 'BinWidth',1)