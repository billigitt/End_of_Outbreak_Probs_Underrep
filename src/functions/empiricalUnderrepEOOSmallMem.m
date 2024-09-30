function output = empiricalUnderrepEOOSmallMem(C_data, RpreERT, RERT, rho, w, ...
    numSamples, tERTArrival, tEndAnalysis, Tinf, experimentsPerSample)

%Function to generate estimate of probability of future cases by simulation. Function inputs generally are understandable in context of paper. w is the serial interval, numSamples is the number of simulated epidemcis (underlying the C_data), and experimentsPerSample is the number of simulations going forward up to Tinf (to observe future cases) for each simulated peidemic. 

%This function is optimised to use less memory.

C_1 = C_data(1);
T = length(C_data);

I_1 = revBin(C_1, rho, numSamples); %generates numSamples X 1 vector of possible values of I_1.

matrixI(1, :) = I_1;

f = waitbar(0,'1','Name','Calculating risks (by simulation)...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

setappdata(f,'canceling',0);

riskWithdrawERT = zeros(tEndAnalysis, 1);

cellOfMatrices = cell(tEndAnalysis, 1);
cellOfMatrices{1} = I_1;


for t = 2:tEndAnalysis

    if getappdata(f,'canceling')
        break
    end

    waitbar((t-2)/(tEndAnalysis-1),f,sprintf("t = "+t+" (of "+2+" to "+tEndAnalysis+")"))


    if (t > (T+1))

        trueC_t = 0;

    else

        trueC_t = C_data(t-1);

    end

    if (t ~= 2)

    numNewIs = 0;
    updatedIs = [];

    R = RpreERT + (RERT - RpreERT)*((t-1) >= tERTArrival); %Decides value of R based on t-1

    while numNewIs < numSamples % could improve this loop so we do a larger batch of simulations on each run

        I_t = renewalEqn(matrixI, w, R); %outputs 1 X numSamples vector with possible values of I_t
        C_t = binornd(I_t, rho*ones(1, numSamples));

        idxAccept = (C_t == trueC_t);

        updatedIs = [updatedIs [matrixI(:, idxAccept); I_t(idxAccept)]];
        numNewIs = numNewIs + sum(idxAccept);

    end
    
    updatedIs = updatedIs(:, 1:numSamples);
    matrixI = updatedIs;%used in the case that we actually get more acceptances than we want
    
    else

        matrixI = I_1';

    end

    cellOfMatrices{t} = matrixI;

    experimentsZeroFutureTrueCases = 0;


    %Memory cheap way to empirically calculate probability of no more cases
    for i = 1:experimentsPerSample

        idxOfInterest = true(1, numSamples);

        for tt = t:Tinf

            updatedIs = renewalEqn(matrixI, w((tt-t+1):end), RpreERT); %NB: Here we just shift the serial interval forward since the first values of w are just mutiplied by 0.
            updatedIs(~idxOfInterest) = 1; %Make previously excluded ones incidence sets non-zero again.

            idxOfInterest(updatedIs ~= 0) = 0;

        end

        experimentsZeroFutureTrueCases = sum(idxOfInterest) + experimentsZeroFutureTrueCases;

    end

    propsOfZeros = experimentsZeroFutureTrueCases/(numSamples*experimentsPerSample);

    riskWithdrawERT(t) = 1 - propsOfZeros;

end

delete(f)


output = struct('riskWithdrawERT', riskWithdrawERT, 'cellOfMatrices', cellOfMatrices);

end
