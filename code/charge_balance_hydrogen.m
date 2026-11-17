function y=charge_balance_hydrogen(x,Ph2,Ph2o,Cd)
global Kph2 Kpo2 Kpeh Kph2o
Cv=x;
%Kph2o=7.09*1e-2;
a=1/Kph2o/Ph2o/Cv;
b=1;
c=Cv-3;
Coh=(-b+sqrt(b^2-4*a*c))/(2*a);
Co=1/Kph2o/Ph2o*Coh^2/Cv;
%Kph2o=4.74*1e-10;
%Co=Kpeh^2/Kph2^2/Kpo2/Cv/sqrt(Po2)/Ph2;
%Kpeh=2.55*1e-14;
Ce=Kph2*sqrt(Ph2)*Co/Coh;
Ch=Kpeh/Ce;
y=2*Cv+Coh+Ch-Ce-Cd;
end