function gamma = gamma_tMostBasic2D(I, w, t)

k = length(w);
T = size(I, 1);
numIncSets = size(I, 2);
gamma = zeros(1, numIncSets);
paddedw = [w; zeros(T-k-1, 1)];

for j = t:T

    gamma = gamma + paddedw(1:(j-1))'*I((j-1):-1:1, :);

end

end