function [CeH,ChH,CohH,CvH,PO2]=boundary_oxygen(ic,TC,Rp,Ph2o,Po2,Cd,num_test)
global R F Kph2o Kpo2 Kpeh Kph2
%constants
R=8.314; %J/(mol*k)
F=96485.0; %C/mol
%molar volume
Vm=5.1323*1e-5; %% m^3/mol
%cI=1800*Vm; %750
%cI=100*Vm;
%%Yosune initial guess
%cI=100*Vm;%600C
cI=100*Vm; %750C
%temperature
TK=TC+273.15;
%chemical reaction constants Robert Kee
% deltaH_O2=-135*10^3; %J/mol
% deltaH_H2O=-93.3*10^3; %J/mol
% deltaH_EH=266.2*10^3; %J/mol
% deltaH_H2=164.35*10^3; %J/mol
% 
% deltaS_O2=-130; %J/mol
% deltaS_H2O=-113.2; %J/mol
% deltaS_EH=0.0; %J/mol
% deltaS_H2=-17.85; %J/mol

%chemical reaction constants Yosune
deltaH_O2=-46*10^3; %J/mol
deltaH_H2O=-62*10^3; %J/mol
deltaH_EH=266.2*10^3; %J/mol
deltaH_H2=-132*10^3; %J/mol

deltaS_O2=-77; %J/mol
deltaS_H2O=-38; %J/mol
deltaS_EH=0.0; %J/mol
deltaS_H2=-8; %J/mol
% deltaH_O2=-36.08*10^3; %J/mol
% deltaH_H2O=-80.3*10^3; %J/mol
% deltaH_EH=266.2*10^3; %J/mol
% deltaH_H2=-146.01*10^3; %J/mol
% 
% deltaS_O2=-125.48; %J/mol
% deltaS_H2O=-71.78; %J/mol
% deltaS_EH=0.0; %J/mol
% deltaS_H2=-0.89; %J/mol

deltaG_O2=deltaH_O2-TK*deltaS_O2;
deltaG_H2O=deltaH_H2O-TK*deltaS_H2O;
deltaG_EH=deltaH_EH-TK*deltaS_EH;
deltaG_H2=deltaH_H2-TK*deltaS_H2;

% Kpo2=exp(-deltaG_O2/R/TK);
% %Kpo2=0.02;
% Kph2o=exp(-deltaG_H2O/R/TK);
%Kpeh=exp(-deltaG_EH/R/TK);
% Kpeh=sqrt(3.73e-7*8.17e-26);
% %Kpeh=8.8487e-10;
% % %Kph2=exp(-deltaG_H2/R/T);
% % 
% % %%Yosune Kph2
% % Kph2=exp(-deltaG_H2/R/TK)*Kpeh;
% Kpeh=100*exp(-deltaG_EH/R/TK)
% % %Kph2=exp(-deltaG_H2/R/T);
% % 
% % %%Yosune Kph2
% % Kph2=exp(-deltaG_H2/R/TK)*Kpeh;
% %deltaGH2O=-245.45*10^3-TK*(-52.51);
% deltaGH2O=-245.45*10^3-TK*(-52.51);
% Kph2o=14.8;
% %Kpeh=sqrt(0.0029*2.7*1e-16);
% Kpo2=2.55e-5;
% Kph2=sqrt(exp(-deltaGH2O/R/TK)*Kph2o*(Kpeh)^2/Kpo2);

% Kpeh=2.06*1e-14;
% Kph2o=1.68;
% Kph2=4.14*1e-8;
% Kpo2=0.3679;
%Kpeh=1.88e-5;
% %Kpeh=sqrt((Kpo2/5.7)*1e-17);
% Kph2=sqrt(exp(-deltaGH2O/R/TK)*Kph2o*(Kpeh)^2/Kpo2);
% % Kpeh=100*exp(-deltaG_EH/R/TK);
% % % %Kph2=exp(-deltaG_H2/R/T);
% % % 
% % % %%Yosune Kph2
% % % Kph2=exp(-deltaG_H2/R/TK)*Kpeh;
% % %deltaGH2O=-245.45*10^3-TK*(-52.51);
% deltaGH2O=-245.45*10^3-TK*(-52.51);
% Kph2o=0.5;
% % %Kpeh=sqrt(0.0029*2.7*1e-16);
%  Kpo2=1.81e-6;
%  Kph2=sqrt(exp(-deltaGH2O/R/TK)*Kph2o*(Kpeh)^2/Kpo2)
%Kph2=7.8e-7;
%Kpo2=1.81e-6;
load('example_parameters_4411.mat');
Kph2o=P(5,num_test);
Kpo2=P(6,num_test);
Kph2=P(7,num_test);
Kpeh=P(8,num_test);

%%%function give initial guess
cI=initial_guess_oxygen(Ph2o,Po2,Cd);
%%%Check initial guess
% Solve the oxygen-electrode boundary charge-balance equation.
% The root is computed with fzero.
%%OCVCV
fun=@(x) charge_balance_oxygen(x,Ph2o,Po2,Cd);
options = optimset('Display','off'); 
[xx fval exitflag output] = fzero(fun,cI,options);
if isnan(xx) == 1
    [xx fval exitflag output]=fzero(fun,2e-3,options);
end
    


Cv=xx;
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

Ph2=(Coh*Ce/(Co*Kph2))^2;
%%EOCV
%EOCV=-R*T/F*log(Ce/Ce0);



CeH=Ce*exp(asinh(ic*Rp*F/R/TK));
ChH=Kpeh/CeH;
%fun=@(x) H2OCH(x,CeH,ChH,Cd,Ph2);
%options = optimset('Display','iter'); 
%[xx fval exitflag output] = fzero(fun,Coh,options);
%CohH=xx
%CvH=0.5*(Cd+CeH-ChH-CohH)
%CoH=3-CvH-CohH;
%Cd=0.1;
CvH=CeH+Cd-ChH-(3-Co);
CohH=3-Co-CvH;
%Kpo2=1.26;
PO2=(Co*ChH^2/(Kpo2*CvH))^2;
end