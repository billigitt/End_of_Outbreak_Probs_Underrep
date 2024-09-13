function output = empiricalUnderrepEOO(C_data, RpreERT, RERT, rho, w, ...
    numSamples, tERTArrival, tEndAnalysis, Tinf, experimentsPerSample)

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


    if (t > T)

        trueC_t = 0;

    else

        trueC_t = C_data(t);

    end

    numNewIs = 0;
    updatedIs = [];

    R = RpreERT + (RERT - RpreERT)*(t >= tERTArrival); %Decides value of R based on t

    while numNewIs < numSamples % could improve this loop so we do a larger batch of simulations on each run

        I_t = renewalEqn(matrixI, w, R); %outputs 1 X numSamples vector with possible values of I_t
        C_t = binornd(I_t, rho*ones(1, numSamples));

        idxAccept = (C_t == trueC_t);

        updatedIs = [updatedIs [matrixI(:, idxAccept); I_t(idxAccept)]];
        numNewIs = numNewIs + sum(idxAccept);

    end

    matrixI = updatedIs(:, 1:numSamples); %used in the case that we actually get more acceptances than we want

    cellOfMatrices{t} = matrixI;

    experimentMatrixI = repmat(matrixI, 1, experimentsPerSample);

    for tt = (t+1):Tinf

        I_tt = renewalEqn(experimentMatrixI, w, RpreERT);
        idx0s = (I_tt == 0);
        experimentMatrixI = [experimentMatrixI(:, idx0s); I_tt(idx0s)]; %rejects simulations that do not give 0s, making computation more efficient

        if (numel(experimentMatrixI) == 0)
            
            experimentsZeroFutureTrueCases = 0;

            break

        end

    end

    experimentsZeroFutureTrueCases = size(experimentMatrixI, 2);

    propsOfZeros = experimentsZeroFutureTrueCases/(numSamples*experimentsPerSample);

    riskWithdrawERT(t) = 1 - propsOfZeros;

end

delete(f)


output = struct('riskWithdrawERT', riskWithdrawERT, 'cellOfMatrices', cellOfMatrices);

end