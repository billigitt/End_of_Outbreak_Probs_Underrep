function gamma = gamma_tFaster(I, w, t)

T = length(I);
k = length(w);

if k < (T-1)

    paddedw = [w; zeros(T-k-1, 1)];

end

matrixOfIs = repmat(I(1:T-1)', T-1, 1);
matrixOfIs = full(spdiags(zeros(T-1, T+t), [1:(T-1) (t-T-1):-1:(1-T)], matrixOfIs));



gamma = sum(matrixOfIs, 2)'*paddedw(1:(T-1));

end