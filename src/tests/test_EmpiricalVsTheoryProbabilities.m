%test empiricalInf0s with theory

I = [1 0 1]';

w = [0.1 0.2 0.3 0.4]';

R = 1;
numExp = 600000;
Tinf = 25;

prob0sEmp = empiricalInf0s(I, w, R, numExp, Tinf)

gamma = gamma_tMostBasic([I; zeros(Tinf - length(I), 1)], w, length(I)+1);
prob0sTheory = exp(-R*gamma)