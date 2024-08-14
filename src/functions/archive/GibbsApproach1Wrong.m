function probMoreCases = GibbsApproach1Wrong(C, w, rho, GibbsSamples, RpreERT, RERT, tERTDeployed)

I = round((1+C)/rho);
T = length(I);

probGibbsSample = zeros(T, GibbsSamples);
gammaTmp = zeros(1, T);

for j = 1:GibbsSamples

    for t = 1:T

        rangeI = C(t):round(2*(1+C(t))/rho);

        pmfTmp = likelihoodOfTrueIncidence(I, C, w, RpreERT, RERT, rangeI, t, tERTDeployed, rho);
        pmfTmp = round(pmfTmp/sum(pmfTmp), 10); %does this still add to 1? This is currently the PMF
        % for the hidden number of infections.

        I(t) = randsample(length(pmfTmp), 1, true, pmfTmp) + C(t) - 1; % adding C(t) - 1 makes it a sample of the true # of infections

    end

    for t = 1:T

        gammaTmp(t) = gamma_tMostBasic(I, w, t);

    end

    probGibbsSample(:, j) = 1 - exp(-RERT*gammaTmp); % This is used as opposed to
    %riskWithdrawERT.m because we are not using a distributional estimate.

end

probMoreCases = mean(probGibbsSample, 2);

end