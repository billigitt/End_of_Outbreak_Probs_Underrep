function likelihood = likelihoodOfTrueIncidenceSlow(I, C, w, RpreERT, RERT, rangeI, t, tERTDeployed, rho)

lengthRangeI = length(rangeI);

T = length(I);

temporalR = [ones(tERTDeployed - 2, 1)*RpreERT; ones(T - tERTDeployed + 1, 1)*RERT];
temporalR = repmat(temporalR, 1, lengthRangeI);

I3D = zeros(T-1, T-1, lengthRangeI); %This matrix is set up to make the multiplication in calculating gamma
% as computationally efficient as possible
I2D = repmat(I(2:T), 1, lengthRangeI);
%I2D does not need I(1) at all because it is only used in the summand part
%of the likelihood, and this runs from k = 2 to T.
if (t ~= 1)

    I2D(t-1, :) = rangeI; %Note we use t-1 here because the matrix index no longer lines up with the time, t.

end

%From gamma equation, we know that the serial interval only ever needs to
%be as long as the incidence time-series length, minus 1 (as a maximum).

if length(w) <= (T-1)

    paddedw = [w; zeros(T - 1 - length(w), 1)];

else

    paddedw = w(1:(T-1));

end

for j = 1:lengthRangeI

    I(t) = rangeI(j);

    I3D(:, :, j) = toeplitz([I(1) zeros(1, T-2)], I(1:T-1)); %rows are direction of dot product with serial
    % cols are the index k in the equation, 3rd dim are for different
    % values of I(t)

end

gamma = squeeze(pagemtimes(paddedw', I3D));

summand = I2D.*log(temporalR.*gamma) + log(1-rho)*I2D - temporalR.*gamma - log(factorial(I2D - C(2:T)));

if (t ~= 1)

%     likelihood = sum(summand, 1) + factorial(I(1)) + I(1)*log(1-rho) + log(factorial(I(1) - C(1)));
    loglikelihood = sum(summand, 1) + I(1)*log(1-rho);

    %not computing the factorial, and instead writing out the sum of the
    %logs makes rounding errors less likely. NB: in doing this, some of the
    %summations are cancelled out.

    if (C(1) > 0)

        for ii = (I(1) - C(1) + 1):I(1)

            loglikelihood = loglikelihood + log(ii);

        end



    end

else

%     likelihood = sum(summand, 1) + factorial(rangeI) + rangeI*log(1-rho) + log(factorial(rangeI - C(1)));
    logFactorialMatrix = log(hankel(1:C(1), rangeI)); %hankel matrix has k (which is used in summing) on the rows and along columns is teh different possible values of I_1
    loglikelihood = sum(summand, 1) + rangeI*log(1-rho) +sum(logFactorialMatrix, 1);

end

likelihood = exp(loglikelihood);

end