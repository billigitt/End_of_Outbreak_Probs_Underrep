function output = bruteForceSimulationEOO(C_data, RpreERT, RERT, rho, w, ...
    numSamples, tERTArrival, tEndAnalysis, Tinf, experimentsPerSample)

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

    if t <= (length(C_data)+1)

        trueC_2tm1 = C_data(2:(t-1));

    else

        trueC_2tm1 = [C_data(2:end); zeros(t-1-length(C_data), 1)];

    end

    incidenceTrue = [];
    numAccepted = 0;

    if (t > 2)

        while (numAccepted < numSamples)

            matrixI = I_1';



            for i = 2:(t-1)

                R = RpreERT + (RERT - RpreERT)*(i >= tERTArrival); %Decides value of R based on t
                matrixI = [matrixI; renewalEqn(matrixI, w, R)];

            end



            matrixCAfter1 = binornd(matrixI(2:end, :), rho);

            logicalMat = (matrixCAfter1 == trueC_2tm1);
            totalEqualVec = sum(logicalMat, 1);
            idxEquate = (totalEqualVec == (t-1));

            if sum(idxEquate) >= 1

                incidenceTrue = [incidenceTrue matrixI(:, idxEquate)];

            end

            numAccepted = numAccepted + sum(idxEquate);

        end

        incidenceTrue = incidenceTrue(:, 1:numSamples);

    else

        incidenceTrue = I_1';

    end

    

    %Memory cheap way to empirically calculate probability of no more cases
    experimentsZeroFutureTrueCases = 0;

    for i = 1:experimentsPerSample

        idxOfInterest = true(1, numSamples);

        for tt = t:Tinf

            updatedIs = renewalEqn(incidenceTrue, w((tt-t+1):end), RpreERT); %NB: Here we just shift the serial interval forward since the first values of w are just mutiplied by 0.
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