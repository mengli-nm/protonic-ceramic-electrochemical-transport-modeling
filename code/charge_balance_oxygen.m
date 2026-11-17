function y=charge_balance_oxygen(x,Ph2o,Po2,Cd)
global R F Kph2o Kpo2 Kpeh
Cv=x;
%Kph2o=7.09*1e-2;
a=1/Kph2o/Ph2o/Cv;
b=1;
c=Cv-3;
Coh=(-b+sqrt(b^2-4*a*c))/(2*a);
%Kph2o=7.09*1e-2;
Co=1/Kph2o/Ph2o*Coh^2/Cv;
%Kpeh=2.55*1e-14;
%Kpo2=1.26;
Ch=sqrt(Cv*Kpo2*sqrt(Po2)/Co);
%Ce=Kpeh*sqrt(Ph2)*Co/Coh;
Ce=Kpeh/Ch;
y=2*Cv+Coh+Ch-Ce-Cd;
end
