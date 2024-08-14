function probMoreCases = GibbsApproach1(C, w, rho, GibbsSamples, RpreERT, RERT, tERTDeployed, tEndCalc)

%Most basic approach: We assume that we know the Rt values exactly, and use
%this to inform our inference.

%tEndCalc is the time point up to which we want to calculate probability
%estimates

f = waitbar(0,'1','Name','Calculating risks...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

setappdata(f,'canceling',0);

I = round((1+C)/rho); %'basic' estimate of I
T = length(I);

probGibbsSample = zeros(T, GibbsSamples);

Tinf = 300; %perhaps make this an input?

for t = tERTDeployed:tEndCalc

    if getappdata(f,'canceling')
        break
    end

    waitbar((t-tERTDeployed)/(tEndCalc-tERTDeployed),f,sprintf("t = "+t+" (of "+tERTDeployed+" to "+tEndCalc+")"))

    tmpI = I(1:t); %Always goes back to the 'basic' estimate of I

    for j = 1:GibbsSamples

        for i = 1:t

            rangeI = C(i):round(2*(1+C(i))/rho); %perhaps make this a larger range?

            pmfTmp = likelihoodOfTrueIncidence(tmpI, C(1:t), w, RpreERT, RERT, rangeI, i, tERTDeployed, rho);
            pmfTmp = round(pmfTmp/sum(pmfTmp), 10); %does this still add to 1? This is currently the PMF
            % for the hidden number of infections.
            cmfTmp = cumsum(pmfTmp);

            %tmpI(i) = randsample(length(pmfTmp), 1, true, pmfTmp) + C(i) - 1; % adding C(t) - 1 makes it a sample of the true # of infections
            tmpI(i) = find(rand <= cmfTmp, 1) + C(i) - 1; %check this line is OK- copied from find part from internet
        end



        gammaTmp_t = gamma_tMostBasic([tmpI; zeros(Tinf-t, 1)], w, t);
        %add on zeros to get probability of no more cases

        probGibbsSample(t, j) = 1 - exp(-RpreERT*gammaTmp_t); % This is used as opposed to
        %riskWithdrawERT.m because we are not using a distributional estimate.

    end

end
delete(f)

probMoreCases = mean(probGibbsSample, 2);

end