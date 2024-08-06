function prob_t = riskWithdrawERT(shape, rate, gamma_t)

prob_t = 1.-(rate./(rate+gamma_t)).^shape;

end