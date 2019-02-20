%Note, this is an equilibrium solution with Temperature conditions at walls
% R2-R1 must be divisible by rstep
T1 = 500;
T2 = 300;
R1 = 0.5;
R2 = 3;

rstep = 0.01;
numR = (R2-R1)/rstep;
r = zeros(numR,1);
Tcurve = r;
r(1) = R1;
Tcurve(1) = T1;

for i = 2:numR
    r(i) = r(i-1) + rstep;
    Tcurve(i) = ((T1-T2)/log(R1/R2))*log(r(i)/R2) + T2;
    
end

plot(r,Tcurve)