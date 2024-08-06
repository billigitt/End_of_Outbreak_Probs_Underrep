function p_t = forceOfInfection(I, w, meanR, t)
%This function calculates the force of infection at time t. it is then used
%(outside of this function) to estimate the relative probability that a
%case has been under-reported.

incidenceConsidered = I((t-1):-1:1);

lengthSI = length(w);
lengthIncidence = length(incidenceConsidered);

%Get them the same length
if (lengthSI >= lengthIncidence)

    w = w(1:lengthIncidence);

else

    incidenceConsidered = incidenceConsidered(1:lengthSI);

end

p_t = meanR*dot(w, incidenceConsidered);

end