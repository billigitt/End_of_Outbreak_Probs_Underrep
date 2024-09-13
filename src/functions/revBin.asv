function I_1 = revBin(C_1, rho, numSamples)

rangeI_1 = C_1:ceil(5*(1 + C_1)/rho);

if C_1 == 0

    prob = (1-rho).^rangeI_1;
    
else

    logProb = hankel((1:C_1)', rangeI_1);
    logProb = log(logProb);
    logProb = sum(logProb, 1) + rangeI_1*log(1-rho);
    prob = exp(logProb);

end

prob = prob/sum(prob);

I_1 = randsample(length(prob), numSamples, true, prob) + C_1 - 1; % adding C(t) - 1;

end