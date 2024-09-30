function prob = empiricalInf0s(I, w, R, numExp, Tinf)

matrixI = repmat(I, 1, numExp);

t = length(I);

for tt = t+1:Tinf

    matrixI = [matrixI; renewalEqn(matrixI, w', R)];

end

totalFutureCases = sum(matrixI(t+1:Tinf, :));

prob = sum(totalFutureCases == 0)/numExp;

end