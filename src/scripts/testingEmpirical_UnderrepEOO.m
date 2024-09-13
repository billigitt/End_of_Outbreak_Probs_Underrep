% generate incidence

addpath('../functions')

w = [0.4 0.3 0.2 0.1];
I_1 = 10;

T = 10;
tERTArrival = 5;

RpreERT = 2;
RERT = 0.5;

I = [I_1; zeros(T-1, 1)];

for t = 2:T

    R_t = RpreERT + (RERT - RpreERT)*(t>=tERTArrival);

    if length(w)<(t-1)

        wTmp = [w zeros(1, t - 1 - length(w))];

    else

        wTmp = w(1:(t-1));
    
    end

    I(t) = poissrnd(R_t*dot(I((t-1):-1:1)', wTmp));

end

rho = 0.4;

C = binornd(I(t), rho*ones(T, 1));

numSamples = 1e4;
tEndAnalysis = T*2;
Tinf = 1000;
experimentsPerSample = 1;

risks = empiricalUnderrepEOO(C, RpreERT, RERT, rho, w, numSamples, tERTArrival, tEndAnalysis, Tinf, experimentsPerSample)