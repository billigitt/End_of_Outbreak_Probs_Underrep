function gamma = gamma_tMostBasic3D(I, w, t)

% I is now 3D, i.e. time X sample X m, and so the output should be a sample
% X m matrix, instead of a scalar.

k = length(w);
T = size(I, 1);
numSamples = size(I, 2); numMs = size(I, 3);
gamma = zeros(numSamples, numMs);
paddedw = [w' zeros(1, T-k-1)];

for j = t:T

    gamma = gamma + squeeze(pagemtimes(paddedw(1:(j-1)), I((j-1):-1:1, :, :)));

end

end