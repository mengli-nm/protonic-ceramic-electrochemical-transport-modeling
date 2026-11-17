function [CeH2,ChH2,CohH2,CvH2,PO22]=boundary_hydrogen(ic,TC,Rp,Ph2,Ph2o,Cd,num_test)
global R F Kph2o Kpo2 Kpeh Kph2
%constants
R=8.314; %J/(mol*k)
F=96485.0; %C/mol
%molar volume
Vm=5.1323*1e-5; %% m^3/mol
%cI=2000*Vm; %750
%cI=50*Vm;

%Yosune initial guess
%cI=10000*Vm; %600C
cI=2000*Vm;
%temperature
TK=TC+273;
%chemical reaction constants
% deltaH_O2=-135*10^3; %J/mol
% deltaH_H2O=-93.3*10^3; %J/mol
% deltaH_EH=266.2*10^3; %J/mol
% deltaH_H2=164.35*10^3; %J/mol
% 
% deltaS_O2=-130; %J/mol
% deltaS_H2O=-113.2; %J/mol
% deltaS_EH=0.0; %J/mol
% deltaS_H2=-17.85; %J/mol

%chemical reaction constants yosune
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
% Kpeh=exp(-deltaG_EH/R/TK);
% % % deltaGH2O=-248.11*10^3-TK*(-55.48);
% % % Kpeh=sqrt(3.73e-7*8.17e-26);
% % % Kph2o=25.89;
% % % Kpo2=3.37e-7;
% % % Kph2=sqrt(exp(-deltaGH2O/R/TK)*Kph2o*(Kpeh)^2/Kpo2);
%Kpeh=8.8487e-10;
% %Kph2=exp(-deltaG_H2/R/T);
% %Yosune Kph2
% Kph2=exp(-deltaG_H2/R/TK)*Kpeh;
% deltaGH2O=-248.11*10^3-TK*(-55.48);
% %deltaGH2O=-245.45*10^3-TK*(-52.51);
% %Kph2o=1/0.23;
% Kph2o=25.89;
% %Kph2o=14.8;
% %Kpeh=sqrt(0.0029*2.7*1e-16);
% %Kpo2=0.0029^2*1/0.23;
% %Kpo2=0.0029^2*Kph2o;
% Kpo2=3.37e-7;
%%Kpeh=1.88e-15;
% %Kpeh=sqrt((Kpo2/5.7)*1e-17);
% Kph2=sqrt(exp(-deltaGH2O/R/TK)*Kph2o*(Kpeh)^2/Kpo2);
% % Kpeh=100*exp(-deltaG_EH/R/TK);
% % % %Kph2=exp(-deltaG_H2/R/T);
% % % 
% % % %%Yosune Kph2
% % % Kph2=exp(-deltaG_H2/R/TK)*Kpeh;
% % %deltaGH2O=-245.45*10^3-TK*(-52.51);
%% deltaGH2O=-245.45*10^3-TK*(-52.51);
%% Kph2o=0.5;
% % %Kpeh=sqrt(0.0029*2.7*1e-16);
%  Kpo2=1.81e-6;
%  Kph2=sqrt(exp(-deltaGH2O/R/TK)*Kph2o*(Kpeh)^2/Kpo2)
load('example_parameters_4411.mat');
Kph2o=P(5,num_test);
Kpo2=P(6,num_test);
Kph2=P(7,num_test);
Kpeh=P(8,num_test);

%%Kph2=7.8e-7;
%%Kpo2=1.81e-6;

%%%Check initial guess
% xx=(1e-3:1e-3:1);
% for ii=1:1:1000
%     yy(ii)=charge_balance_hydrogen(xx(ii),Ph2,Ph2o,Cd);
% end
% plot(xx,yy)
%%%function give initial guess
cI=initial_guess_hydrogen(Ph2,Ph2o,Cd);

fun=@(x) charge_balance_hydrogen(x,Ph2,Ph2o,Cd);
options = optimset('Display','off'); 
[xx fval exitflag output] = fzero(fun,cI,options);
Cv=xx;
a=1/Kph2o/Ph2o/Cv;
b=1;
c=Cv-3;
Coh=(-b+sqrt(b^2-4*a*c))/(2*a);
Co=1/Kph2o/Ph2o*Coh^2/Cv;
%Kpeh=2.55*1e-14;
Ce=Kph2*sqrt(Ph2)*Co/Coh;
Ch=Kpeh/Ce;

%after charging
CeH2=Ce*exp(-asinh(ic*Rp*F/R/TK));
%size(CeH2)
%size(Kpeh)
ChH2=Kpeh/CeH2;
CvH2=CeH2+Cd-ChH2-(3-Co);
CohH2=3-Co-CvH2;
PO22=(Co*ChH2^2/(Kpo2*CvH2))^2
end