function I = renewalEqn(matrixI, w, R)

T = size(matrixI, 1);
lengthSI = length(w);

if (lengthSI < T)

    w = [w zeros(1, T - lengthSI)];

else

    w = w(1:T);

end

gammaVec = fliplr(w)*matrixI;
infPot = gammaVec*R;

I = poissrnd(infPot);

end