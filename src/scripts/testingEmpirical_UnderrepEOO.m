% generate incidence

% Add functions to filepath
addpath('../functions')

% Speccify the serial interval
w = [0.4 0.3 0.2 0.1];

% Set number of cases in initial timestep
I_1 = 10;

% Set total number of timesteps
T = 10;

% Set timestep that Emmergency Response Team (ERT) arrives
tERTArrival = 5;

% Set R values during two time periods
RpreERT = 2; % Prior to ERT time being active
RERT = 0.5; % During ERT active period

% Initialise vector for incidence count per timestep
I = [I_1; zeros(T-1, 1)];

% Loop over each timestep and populate incidence count for that timestep
for t = 2:T
    
    % Get R_t value based on whether in preERT or ERT arrival period
    R_t = RpreERT + (RERT - RpreERT)*(t>=tERTArrival);
    
    % Construct temporary variable that contains serial interval vector and
    % is of length t-1 (pad serial interval "w" with zeros if needed)
    if length(w)<(t-1)

        wTmp = [w zeros(1, t - 1 - length(w))];

    else

        wTmp = w(1:(t-1));
    
    end

    % Generate incidence count for timestep t
    I(t) = poissrnd(R_t*dot(I((t-1):-1:1)', wTmp));

end

% Specify reporting factor rho
% Generate reported case count for each timestep (binomial sampling)
rho = 0.4;
C = binornd(I(t), rho*ones(T, 1));

% Using the reported cases C and reporting rate rho, numerically estimate
% the end of outbreak probability
numSamples = 1e4;
tEndAnalysis = T*2;
Tinf = 1000;
experimentsPerSample = 1;
risks = empiricalUnderrepEOO(C, RpreERT, RERT, rho, w, numSamples, tERTArrival, tEndAnalysis, Tinf, experimentsPerSample)