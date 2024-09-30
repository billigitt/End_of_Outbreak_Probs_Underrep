function [shape, rate] = posteriorRtDuringERT(I, w, tDeploy, tWithdraw)

%function to find posterior distribution of Rt over entire time series
%between Ebola Response Team (ERT) being deployed and being withdrawn.


shape = 1 + sum(I(tDeploy:tWithdraw));

rate = 0;
for t = (tDeploy:tWithdraw) 

    rate = rate + forceOfInfection(I, w, 1, t);

end

end