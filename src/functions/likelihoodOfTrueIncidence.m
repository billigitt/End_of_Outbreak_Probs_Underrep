function likelihood = likelihoodOfTrueIncidence(I, C, w, RpreERT, RERT, rangeI, t, tERTDeployed, rho)

lengthRangeI = length(rangeI);

T = length(I);

temporalR = [ones(tERTDeployed - 2, 1)*RpreERT; ones(T - tERTDeployed + 1, 1)*RERT];
temporalR = repmat(temporalR, 1, lengthRangeI);

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

firstRow = [I(1); zeros(T-2, 1)];

%perturbationMatrix is just how the contribution to gamma changes when
%different values of I(t) are used (starting from I(t) = C(t) up to the
%max considered).
%firstRow and firstColumn are used to generated the Toeplitz matrix, which
%enables computationally efficient matrix-vector multiplication to
%calculate gamma when I(t) = C(t). We then use the pertrbation matrix to
%simply add on all the xtra contributions for different values of I(t).
if t == 1

    firstColumn = [C(1); I(2:(T-1))];
    firstRow(1) = C(1);
    perturbationMatrix = paddedw.*(0:(lengthRangeI-1));

elseif t == (T-1)

    firstColumn = [I(1:(T-2)); C(T-1)];
    perturbationMatrix = [zeros(T-2, lengthRangeI); paddedw(1).*(0:(lengthRangeI-1))];

elseif t == T

    firstColumn = I(1:(T-1));
    perturbationMatrix = zeros(T-1, lengthRangeI);

else

    firstColumn = [I(1:(t-1)); C(t); I((t+1):(T-1))];
    perturbationMatrix = [zeros(t-1, lengthRangeI); paddedw(1:(T-t)).*(0:(lengthRangeI-1))];


end

gamma = toeplitz(firstColumn, firstRow)*paddedw;
gamma = repmat(gamma, 1, lengthRangeI);
gamma = gamma + perturbationMatrix;

summand = I2D.*log(temporalR.*gamma);
summand = summand + log(1-rho)*I2D;
summand = summand - temporalR.*gamma; 
summand = summand - log(factorial(I2D - C(2:T))); %this line could be sped-up if we expand logarithm part?

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