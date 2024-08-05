function gamma = gamma_tMostBasic(I, w, t)

k = length(w);
T = length(I);
gamma = 0;
paddedw = [w; zeros(T-k-1, 1)];

for j = t:T

    gamma = gamma + dot(I((j-1):-1:1), paddedw(1:(j-1)));

end

end