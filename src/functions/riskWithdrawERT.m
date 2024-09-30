function prob_t = riskWithdrawERT(shape, rate, gamma_t)

%Taken from Thompson et al (Nature communications, 2024) to calculate the probability of future cases from a distributional (gamma) estimate of the reproduction number.

prob_t = 1.-(rate./(rate+gamma_t)).^shape;

end
