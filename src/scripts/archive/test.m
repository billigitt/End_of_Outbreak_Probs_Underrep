SI_mean_daily = 15.3;
SI_sd_daily = 9.3;

SI_scale_daily = SI_sd_daily^2/SI_mean_daily;
SI_shape_daily = SI_mean_daily/SI_scale_daily;


wDaily = zeros(100,1);
for k = 1:100
    intVals = [k-1:(1/10000):k+1];
    funcVals = (1 - abs(intVals - k)).*gampdf(intVals,SI_shape_daily,SI_scale_daily);
    wDaily(k) = trapz(intVals, funcVals);
end

wDaily(1) = wDaily(1) + (1-sum(wDaily));

% Weekly serial interval calc

P = 7;
SI_mean_weekly = 15.3/P;
SI_sd_weekly = 9.3/P;

SI_scale_weekly = SI_sd_weekly^2/SI_mean_weekly;
SI_shape_weekly = SI_mean_weekly/SI_scale_weekly;

wDaily2 = zeros(100,1);
for k = 1:100
    intVals = linspace((k-1)/P, (k+1)/P, 1000);
    funcVals = (1 - P*abs(intVals - (k/P))).*gampdf(intVals,SI_shape_weekly,SI_scale_weekly);
    wDaily2(k) = trapz(intVals, funcVals);
end

wDaily2(1) = wDaily2(1) + (1-sum(wDaily2));