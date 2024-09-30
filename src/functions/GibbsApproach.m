function output = GibbsApproach(C, w, rho, GibbsSamples, RpreERT, RERT, ...
    tERTDeployed, tEndCalc, burnin, thinning)

%Function to generate estimate of probability of future cases by Gibbs Sampling. Function inputs generally are understandable in context of paper. w is the serial interval, GibbsSamples is the number of iterations in the Gibbs sampling method.

%Most basic approach: We assume that we know the Rt values exactly, and use
%this to inform our inference.

%tEndCalc is the time point up to which we want to calculate probability
%estimates

f = waitbar(0,'1','Name','Calculating risks...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

setappdata(f,'canceling',0);

I = round((C)/rho); %'basic' estimate of I
T = length(I);

probGibbsSample = zeros(T, GibbsSamples);

Tinf = 200;

incidenceStore = cell(tEndCalc - tERTDeployed + 1, 1);
pGeweke = zeros(tEndCalc, 1);
pGewekeBuiltIn = zeros(tEndCalc, 1);
zGeweke = zeros(tEndCalc, 1);

for t = tERTDeployed:tEndCalc

    if getappdata(f,'canceling')

        break
        
    end

    waitbar((t-tERTDeployed)/(tEndCalc-tERTDeployed),f,sprintf("t = "+t+" (of "+tERTDeployed+" to "+tEndCalc+")"))

    tmpI = I(1:(t-1)); %Always goes back to the 'basic' estimate of I.

    incidenceStore{t - tERTDeployed + 1} = zeros((t-1), GibbsSamples);

    for j = 1:GibbsSamples

        for i = 1:(t-1)

            rangeI = C(i):ceil(5*(1+C(i))/rho); 

            if (rho == 1)

                tmpI(i) = C(i);

            else

                pmfTmp = likelihoodOfTrueIncidence(tmpI, C(1:(t-1)), w, RpreERT, RERT, rangeI, i, tERTDeployed, rho);
                pmfTmp = round(pmfTmp/sum(pmfTmp), 10); 
                % for the hidden number of infections.
                cmfTmp = cumsum(pmfTmp);

                tmpI(i) = find(rand <= cmfTmp, 1) + C(i) - 1; 

            end

        end

        incidenceStore{t - tERTDeployed + 1}(:, j) = tmpI;

        gammaTmp_t = gamma_tMostBasic([tmpI; zeros(Tinf-t+1, 1)], w, t);
        %add on zeros to get probability of no more cases

        probGibbsSample(t, j) = 1 - exp(-RpreERT*gammaTmp_t); % This is used as opposed to
        %riskWithdrawERT.m because we are not using a distributional estimate.

    end

    [zGeweke(t), pGewekeBuiltIn(t)] = geweke(probGibbsSample(t, burnin + 1:thinning:end)');
    [~, pGeweke(t) ] = ztest(zGeweke(t), 0, 1);
end
delete(f)

probMoreCases = mean(probGibbsSample(:, burnin + 1:thinning:end), 2);


output = struct('probMoreCases', probMoreCases, 'incidenceStore', incidenceStore, 'pGeweke', pGeweke, 'pGewekeBuiltIn', pGewekeBuiltIn, 'probGibbsSample', probGibbsSample);

end
