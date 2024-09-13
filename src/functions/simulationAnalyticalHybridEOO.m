function output = simulationAnalyticalHybridEOO(C_data, RpreERT, RERT, rho, w, ...
    numSamples, tERTArrival, tEndAnalysis, Tinf)

C_1 = C_data(1);
T = length(C_data);

I_1 = revBin(C_1, rho, numSamples); %generates numSamples X 1 vector of possible values of I_1.

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

    if t <= length(C_data)

        trueC_2t = C_data(2:t);

    else

        trueC_2t = [C_data(2:end); zeros(t-length(C_data), 1)];

    end

    incidenceTrue = [];
    numAccepted = 0;

    while (numAccepted < numSamples)

        matrixI = I_1';

        for i = 2:t

            R = RpreERT + (RERT - RpreERT)*(i >= tERTArrival); %Decides value of R based on t
            matrixI = [matrixI; renewalEqn(matrixI, w, R)];

        end

        matrixCAfter1 = binornd(matrixI(2:end, :), rho);

        logicalMat = (matrixCAfter1 == trueC_2t);
        totalEqualVec = sum(logicalMat, 1);
        idxEquate = (totalEqualVec == (t-1));

        if sum(idxEquate) >= 1

            incidenceTrue = [incidenceTrue matrixI(:, idxEquate)];

        end

        numAccepted = numAccepted + sum(idxEquate);

    end

    incidenceTrue = incidenceTrue(:, 1:numSamples);
    incidenceTrue = [incidenceTrue; zeros(Tinf-t, numSamples)];

    gammaTmp = gamma_tMostBasic2D(incidenceTrue, w', t);
    
    riskWithdrawERTAll = 1 - exp(-RpreERT*gammaTmp);

    riskWithdrawERT(t) = mean(riskWithdrawERTAll);

end

delete(f)

output = struct('riskWithdrawERT', riskWithdrawERT, 'cellOfMatrices', cellOfMatrices);

end