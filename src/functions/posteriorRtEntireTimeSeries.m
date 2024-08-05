function [shape, rate] = posteriorRtEntireTimeSeries(IbeforeERT, wCDF)

%function to find posterior distribution of Rt over entire time series
%prior to Ebola Response Team (ERT) being deployed.

%tests: check wCDF sums to 1. check IbeforeERT matches what we have before

T = length(IbeforeERT);
lengthSI = length(wCDF);

shape = 1 + sum(IbeforeERT(2:end));

timesConsidered = min(lengthSI, T-1);
incidenceConsidered = IbeforeERT((T-1):-1:(T-timesConsidered));
wCDFConsidered = wCDF(1:timesConsidered);

rate = sum(incidenceConsidered.*wCDFConsidered);

end