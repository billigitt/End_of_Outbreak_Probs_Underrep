function [shape, rate] = posteriorRtEntireTimeSeries3D(IbeforeERT3D, wCDF)

%function to find posterior distribution of Rt over entire time series
%prior to Ebola Response Team (ERT) being deployed. the incidence input is
%now a 3d input.

%tests: check wCDF sums to 1. check IbeforeERT matches what we have before

T = size(IbeforeERT3D, 1);
lengthSI = length(wCDF);

shape = 1 + squeeze(sum(IbeforeERT3D(2:end, :, :), 1));

timesConsidered = min(lengthSI, T-1);
incidenceConsidered = IbeforeERT3D((T-1):-1:(T-timesConsidered), :, :);
wCDFConsidered = wCDF(1:timesConsidered);

rate = squeeze(pagemtimes(wCDFConsidered', incidenceConsidered));

end