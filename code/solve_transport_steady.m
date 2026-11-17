%%%Solve the charge species distribution in the electrolyte
%function [ieh,Ce,Ch,Coh,Cv,
function [yinit,CeS,ChS,CpS,CvS,CoS,PhiS,CoSOCV,PO2SOCV,PO2S,PH2S,PH2OS,rleak,rleakS,rv,rp,re,rh,icS,V]=solve_transport_steady(icu,Temp,Thick,RpH,RpH2O,Phh2o,Phh2,Psh2o,Pso2,Cd,icguess,num_test)
global ChH CeH CpH CvH ChOH CeOH CpOH CvOH ic T F R NumGrid Vm regupara
global Dh De Dv Dp zh ze zv zp dx epsilon0 epsilonr syp Kpo2 Kpeh Kph2 Kph2o

zh=1.0;
ze=-1.0;
zv=2.0;
zp=1.0;
T=Temp+273.15;

NumGrid=200;
R=8.314; %J/(mol*k)
F=96485.0; %C/mol
%Vm=4.55*1e-5;
Vm=5.1323*1e-5; %%m^3/mol
Tscale=600; %%s
dx=1/(NumGrid); %0.01L

regupara=0.0;

epsilon0=100.0;
%epsilon0=0.0;
%epsilonr=1e-5;
epsilonr=100.0*Vm*R*T/F/F/Thick/Thick;

% Dhd=1.05*1e-6*exp(-97*1e+3/R/T); %m^2/s
% Ded=1.07*1e-7*exp(-90*1e+3/R/T); %m^2/s
% Dvd=1.90*1e-9*exp(-70*1e+3/R/T);
%Dpd=1.09*1e-7*exp(-43.14*1e+3/R/T);
%Dpd=9.9e-11;
%Dpd=4.23e-11;
% Dpd=2.02e-10;
% %Dvd=2.6e-11;
% %Dvd=1.87*1e-7*exp(-60*1e+3/R/T);
% %Dvd=1.93e-11;
% Dvd=5.86e-11;
% 
% %Dhd=1.52e-8;
% %Dhd=3.78e-8;
% Dhd=3.14e-7;
% 
% %Ded=1.07*1e-7*exp(-90*1e+3/R/T);
% Ded=4.36e-7;
% %Ded=9.81e-8;
% %Ded=2.85e-7;

load('example_parameters_4411.mat');


Dpd=P(1,num_test);
Dvd=P(2,num_test);
Dhd=P(3,num_test);
Ded=P(4,num_test);

%%%Scale diffusivities with Dp
Dh(1:NumGrid)=Dhd./Dpd; %dimensionless
De(1:NumGrid)=Ded./Dpd; %dimensionless
Dv(1:NumGrid)=Dvd./Dpd; 
Dp(1:NumGrid)=Dpd./Dpd;

%%%Scale conrent density with 

ic=icu*Thick*Vm/Dpd/F*1e+4; %icu in the dimension of A/cm^2
syp=Thick^2/Dpd/Tscale;
icg=icguess*Thick*Vm/Dpd/F*1e+4; %icuguess in the dimension of A/cm^2

[CeH2,ChH2,CpH2,CvH2,PO2H]=boundary_hydrogen(icguess,Temp,RpH,Phh2,Phh2o,Cd,num_test);
[CeOH2,ChOH2,CpOH2,CvOH2,PO2H2O]=boundary_oxygen(icguess,Temp,RpH2O,Psh2o,Pso2,Cd,num_test);

%%%Scale concentrations with 1/Vm
CeH=sqrt(CeH2); %dimensionless CeH=Ce*Vm
ChH=sqrt(ChH2);
CpH=sqrt(CpH2);
CvH=sqrt(CvH2);

CeOH=sqrt(CeOH2);
ChOH=sqrt(ChOH2);
CpOH=sqrt(CpOH2);
CvOH=sqrt(CvOH2);

PhiH=0.0;

%%%Scale conrent density with 
% Vm=4.55*1e-5; %%m^3/mol
% ic=icu*Thick*Vm/(1.55*1e-8*exp(-43*1e+3/R/T))/F*1e+4; %icu in the dimension of A/cm^2
% 
% syp=Thick^2/(1.55*1e-8*exp(-43*1e+3/R/T))/60;
% Use an inconsistent initial condition to test initialization.
%y0 = [1; 0; 1e-3];
%yp0 = [0; 0;    0];

%find the initial condition with OCV case

%tspan=(0:6:360*5*10*10);
%tspan=[0 360*5*10*10];

% M = [1 0 0
%    0 1 0
%    0 0 0];

% Use the LSODI example tolerances and set the value of df/dy'

%x0(1:5,1:NumGrid)=1.0;
%x0(1,1:NumGrid)=CeH+(0:1:NumGrid-1).*(CeOH-CeH)/(NumGrid-1);
x0(1,1:NumGrid)=ChH+(0:1:NumGrid-1).*(ChOH-ChH)/(NumGrid-1);
x0(2,1:NumGrid)=CvH+(0:1:NumGrid-1).*(CvOH-CvH)/(NumGrid-1);
x0(3,1:NumGrid)=CpH+(0:1:NumGrid-1).*(CpOH-CpH)/(NumGrid-1);
x0(4,1:NumGrid)=0.0;

% semilogy(x0(1,1:NumGrid))
% pause

%%%solve for OCV case to get the initial condition for the transient
%%%equation

fun=@(y) finitial(y,icg);
yinit=fsolve(fun,x0);

yinit=interpolate_initial_state(yinit);

% for ww=2:1:5
%     ww
% [CeH2,ChH2,CpH2,CvH2,PO2H]=boundary_hydrogen(icguess,Temp,RpH,Phh2,10^(-ww),Cd,num_test);
% [CeOH2,ChOH2,CpOH2,CvOH2,PO2H2O]=boundary_oxygen(icguess,Temp,RpH2O,Psh2o,Pso2,Cd,num_test);
% 
% %%%Scale concentrations with 1/Vm
% CeH=sqrt(CeH2); %dimensionless CeH=Ce*Vm
% ChH=sqrt(ChH2);
% CpH=sqrt(CpH2);
% CvH=sqrt(CvH2);
% 
% CeOH=sqrt(CeOH2);
% ChOH=sqrt(ChOH2);
% CpOH=sqrt(CpOH2);
% CvOH=sqrt(CvOH2);
% x0=yinit;
% 
% fun=@(y) finitial(y,icg);
% yinit=fsolve(fun,x0);
% 
% yinit=interpolate_initial_state(yinit);
% end
icguess1=0;
dc=(icguess1-icu)/5;
for zz=1:1:5
    zz
    icguess1=icguess1-dc
    icg=icguess1*Thick*Vm/Dpd/F*1e+4
[CeH2,ChH2,CpH2,CvH2,PO2H]=boundary_hydrogen(icguess1,Temp,RpH,Phh2,Phh2o,Cd,num_test);
[CeOH2,ChOH2,CpOH2,CvOH2,PO2H2O]=boundary_oxygen(icguess1,Temp,RpH2O,Psh2o,Pso2,Cd,num_test);

%%%Scale concentrations with 1/Vm
CeH=sqrt(CeH2); %dimensionless CeH=Ce*Vm
ChH=sqrt(ChH2);
CpH=sqrt(CpH2);
CvH=sqrt(CvH2);

CeOH=sqrt(CeOH2);
ChOH=sqrt(ChOH2);
CpOH=sqrt(CpOH2);
CvOH=sqrt(CvOH2);
x0=yinit;

fun=@(y) finitial(y,icg);
yinit=fsolve(fun,x0);

yinit=interpolate_initial_state(yinit);
end

CoSOCV=3-yinit(2,:).^2-yinit(3,:).^2;
PO2SOCV=(CoSOCV.*((yinit(1,:).^2)).^2./(yinit(2,:).^2)./Kpo2).^2;


% For this problem you cannot fix more than two components of the guesses.
% Fixing the first two components of y0 leads to the same consistent initial
% conditions found by ODE15S in HB1DAE.
%y0(1:5*NumGrid)=0.0;


 
 


[CeH2,ChH2,CpH2,CvH2,PO2H]=boundary_hydrogen(icu,Temp,RpH,Phh2,Phh2o,Cd,num_test);
[CeOH2,ChOH2,CpOH2,CvOH2,PO2H2O]=boundary_oxygen(icu,Temp,RpH2O,Psh2o,Pso2,Cd,num_test);

CeH=(CeH2); %dimensionless CeH=Ce*Vm
ChH=(ChH2);
CpH=(CpH2);
CvH=(CvH2);

CeOH=(CeOH2);
ChOH=(ChOH2);
CpOH=(CpOH2);
CvOH=(CvOH2); 



%y=reformulateinitial(y);
%sz=size(y);
CeS(1:NumGrid+2)=0.0;
ChS(1:NumGrid+2)=0.0;
CvS(1:NumGrid+2)=0.0;
CpS(1:NumGrid+2)=0.0;
PhiS(1:NumGrid+2)=0.0;

ChS(2:NumGrid+1) = yinit(1,1:NumGrid).^2;
CeS(2:NumGrid+1) = Kpeh./(yinit(1,1:NumGrid).^2);
CvS(2:NumGrid+1) = yinit(2,1:NumGrid).^2;
CpS(2:NumGrid+1) = yinit(3,1:NumGrid).^2;
PhiS(2:NumGrid+1) = yinit(4,1:NumGrid);


ChS(1)=ChH;
CeS(1)=CeH;
CvS(1)=CvH;
CpS(1)=CpH;
%PhiS(1:sz(1),1)=PhiH;

ChS(NumGrid+2)=ChOH;
CeS(NumGrid+2)=CeOH;
CvS(NumGrid+2)=CvOH;
CpS(NumGrid+2)=CpOH;
%PhiS(1:sz(1),1)=PhiH;

CoS=3-CvS-CpS;
PO2S=(CoS.*ChS.^2./CvS./Kpo2).^2;
PH2S=(CpS.*CeS./CoS./Kph2).^2;
PH2OS=CpS.^2./CoS./CvS./Kph2o;

xgrid=(2:1:NumGrid-1);
ih(1:NumGrid)=0.0;
ie(1:NumGrid)=0.0;
iv(1:NumGrid)=0.0;
ip(1:NumGrid)=0.0;

ih(xgrid)=-zh*(Dh(xgrid).*(ChS(xgrid+1)-ChS(xgrid-1)))./(2*dx)...
    -zh^2*Dh(xgrid).*ChS(xgrid).*(PhiS(xgrid+1)-PhiS(xgrid-1))./(2*dx);

ie(xgrid)=-ze*(De(xgrid).*(CeS(xgrid+1)-CeS(xgrid-1)))./(2*dx)...
    -ze^2*De(xgrid).*CeS(xgrid).*(PhiS(xgrid+1)-PhiS(xgrid-1))./(2*dx);

iv(xgrid)=-zv*(Dv(xgrid).*(CvS(xgrid+1)-CvS(xgrid-1)))./(2*dx)...
    -zv^2*Dv(xgrid).*CvS(xgrid).*(PhiS(xgrid+1)-PhiS(xgrid-1))./(2*dx);

ip(xgrid)=-zp*(Dp(xgrid).*(CpS(xgrid+1)-CpS(xgrid-1)))./(2*dx)...
    -zp^2*Dp(xgrid).*CpS(xgrid).*(PhiS(xgrid+1)-PhiS(xgrid-1))./(2*dx);

rleak=1./(1+((ie+ih)./(iv+ip)).^(-1));
icS=(ie+ih+iv+ip).*Dpd.*(1/Vm)*F/Thick;
rleakS=(ie+ih)./(ie+ih+iv+ip);
re=ie./(ie+ih+iv+ip);
rh=ih./(ie+ih+iv+ip);
rv=iv./(ie+ih+iv+ip);
rp=ip./(ie+ih+iv+ip);
ChS=ChS./Vm;
CeS=CeS./Vm;
CvS=CvS./Vm;
CpS=CpS./Vm;
CoS=CoS./Vm;
PhiS=PhiS.*R*T/F;

PhiL=PhiS(1);
PhiR=PhiS(201);
CeL=CeS(1);
CeR=CeS(202);

V=(+R*T/F*log(CeL)-R*T/F*log(CeR))+PhiR-PhiL;

filename=['transport_i' num2str(abs(10*icu)) '_H2O' num2str(Psh2o*100) '_BZCYYb4411_steady.mat'];
save(filename,'CeS','ChS','CpS','CvS','CoS','PhiS','PO2S','PH2S','PH2OS','rleakS','rv','rp','re','rh','V');
%Ch=y(1,:);
%Ce=y(2,:);
%Cv=y(3,:);
%Coh=y(4,:);
%phi=y(5,:);

% y(:,2) = 1e4*y(:,2);
% 
% figure;
% semilogx(t,y);
% ylabel('1e4 * y(:,2)');
% title('Robertson DAE problem with a Conservation Law, solved by ODE15I');
% xlabel('This is equivalent to the stiff ODEs coded in HB1ODE.');
end
% --------------------------------------------------------------------------

function resc = ff(t,yc,ypc)
global ChH CeH CpH CvH ChOH CeOH CpOH CvOH  F ic dx NumGrid ...
       Dh De Dp Dv epsilon0 epsilonr zh ze zp zv syp Kpeh Vm regupara
   
y(1:4,1:NumGrid)=0.0;
yp(1:4,1:NumGrid)=0.0;

for ii=1:NumGrid
    for jj=1:4
        y(jj,ii)=yc(jj+4*(ii-1));
        yp(jj,ii)=ypc(jj+4*(ii-1));
    end
end

dxi2=1/(2*dx*dx);
yy(1,1:NumGrid)=Kpeh./(y(1,1:NumGrid)+regupara);
%%grid points in space
xgrid=(2:1:NumGrid-1);
res(1,xgrid)=-((Dh(xgrid+1)+Dh(xgrid)).*(y(1,xgrid+1)-y(1,xgrid))...
    -(Dh(xgrid-1)+Dh(xgrid)).*(y(1,xgrid)-y(1,xgrid-1)))*dxi2...
    -zh*((Dh(xgrid+1).*y(1,xgrid+1)+Dh(xgrid).*y(1,xgrid)).*(y(4,xgrid+1)-y(4,xgrid))...
    -(Dh(xgrid-1).*y(1,xgrid-1)+Dh(xgrid).*y(1,xgrid)).*(y(4,xgrid)-y(4,xgrid-1)))*dxi2;
res(1,xgrid)=zh*res(1,xgrid)+ze.*(-((De(xgrid+1)+De(xgrid)).*(yy(1,xgrid+1)-yy(1,xgrid))...
    -(De(xgrid-1)+De(xgrid)).*(yy(1,xgrid)-yy(1,xgrid-1)))*dxi2...
    -ze*((De(xgrid+1).*yy(1,xgrid+1)+De(xgrid).*yy(1,xgrid)).*(y(4,xgrid+1)-y(4,xgrid))...
    -(De(xgrid-1).*yy(1,xgrid-1)+De(xgrid).*yy(1,xgrid)).*(y(4,xgrid)-y(4,xgrid-1)))*dxi2);
%res(1,xgrid)=syp*(zh.*y(1,xgrid).^2-ze).*yp(1,xgrid)+res(1,xgrid).*y(1,xgrid).^2;
% res(1,xgrid)=syp*(zh-ze*Kpeh./y(1,xgrid).^2).*yp(1,xgrid)...
%     -zh.*((Dh(xgrid+1)+Dh(xgrid)).*(y(1,xgrid+1)-y(1,xgrid))...
%     -(Dh(xgrid-1)+Dh(xgrid)).*(y(1,xgrid)-y(1,xgrid-1)))*dxi2...
%     -zh^2.*((Dh(xgrid+1).*y(1,xgrid+1)+Dh(xgrid).*y(1,xgrid)).*(y(4,xgrid+1)-y(4,xgrid))...
%     -(Dh(xgrid-1).*y(1,xgrid-1)+Dh(xgrid).*y(1,xgrid)).*(y(4,xgrid)-y(4,xgrid-1)))*dxi2...
%     +ze*Kpeh*((De(xgrid+1)./(y(1,xgrid+1).^2)+De(xgrid)./(y(1,xgrid).^2)).*(y(1,xgrid+1)...
%     -y(1,xgrid))-(De(xgrid-1)./(y(1,xgrid-1).^2)+De(xgrid)./(y(1,xgrid).^2)).*(y(1,xgrid)...
%     -y(1,xgrid-1)))*dxi2-ze^2.*Kpeh*((De(xgrid+1)./y(1,xgrid+1)...
%     +De(xgrid)./y(1,xgrid)).*(y(4,xgrid+1)-y(4,xgrid))...
%     -(De(xgrid-1)./y(1,xgrid-1)+De(xgrid)./y(1,xgrid)).*(y(4,xgrid)-y(4,xgrid-1)))*dxi2;

% res(1,xgrid)=syp*yp(1,xgrid)-((Dh(xgrid+1)+Dh(xgrid)).*(y(1,xgrid+1)-y(1,xgrid))...
%     -(Dh(xgrid-1)+Dh(xgrid)).*(y(1,xgrid)-y(1,xgrid-1)))*dxi2...
%     -zh*((Dh(xgrid+1).*y(1,xgrid+1)+Dh(xgrid).*y(1,xgrid)).*(y(5,xgrid+1)-y(5,xgrid))...
%     -(Dh(xgrid-1).*y(1,xgrid-1)+Dh(xgrid).*y(1,xgrid)).*(y(5,xgrid)-y(5,xgrid-1)))*dxi2;
% 
% res(2,xgrid)=syp*yp(2,xgrid)-((De(xgrid+1)+De(xgrid)).*(y(2,xgrid+1)-y(2,xgrid))...
%     -(De(xgrid-1)+De(xgrid)).*(y(2,xgrid)-y(2,xgrid-1)))*dxi2...
%     -ze*((De(xgrid+1).*y(2,xgrid+1)+De(xgrid).*y(2,xgrid)).*(y(5,xgrid+1)-y(5,xgrid))...
%     -(De(xgrid-1).*y(2,xgrid-1)+De(xgrid).*y(2,xgrid)).*(y(5,xgrid)-y(5,xgrid-1)))*dxi2;


res(2,xgrid)=syp*yp(2,xgrid)-((Dv(xgrid+1)+Dv(xgrid)).*(y(2,xgrid+1)-y(2,xgrid))...
    -(Dv(xgrid-1)+Dv(xgrid)).*(y(2,xgrid)-y(2,xgrid-1)))*dxi2...
    -zv*((Dv(xgrid+1).*y(2,xgrid+1)+Dv(xgrid).*y(2,xgrid)).*(y(4,xgrid+1)-y(4,xgrid))...
    -(Dv(xgrid-1).*y(2,xgrid-1)+Dv(xgrid).*y(2,xgrid)).*(y(4,xgrid)-y(4,xgrid-1)))*dxi2;

res(3,xgrid)=syp*yp(3,xgrid)-((Dp(xgrid+1)+Dp(xgrid)).*(y(3,xgrid+1)-y(3,xgrid))...
    -(Dp(xgrid-1)+Dp(xgrid)).*(y(3,xgrid)-y(3,xgrid-1)))*dxi2...
    -zp*((Dp(xgrid+1).*y(3,xgrid+1)+Dp(xgrid).*y(3,xgrid)).*(y(4,xgrid+1)-y(4,xgrid))...
    -(Dp(xgrid-1).*y(3,xgrid-1)+Dp(xgrid).*y(3,xgrid)).*(y(4,xgrid)-y(4,xgrid-1)))*dxi2;

res(4,xgrid)=epsilon0*epsilonr.*(y(4,xgrid+1)-2*y(4,xgrid)+y(4,xgrid-1)).*dxi2*2 ...
     +(zh*y(1,xgrid)+ze*yy(1,xgrid)+zv*y(2,xgrid)+zp*y(3,xgrid)-0.2);


%%%Set boundary condition for phi

 
%%%%boundary conditions
ii=1;
%  sumCF=Dh(ii)*zh*(y(1,ii)-ChH)+De(ii)*ze*(Kpeh/y(1,ii)-CeH) ...
%      +Dv(ii)*zv*(y(2,ii)-CvH)+Dp(ii)*zp*(y(3,ii)-CpH);
%  PhiH=y(4,ii)+2*dx*(ic+sumCF/2/dx)/(Dh(ii)*zh^2*y(1,ii)+De(ii)*ze^2*Kpeh/y(1,ii) ...
%      +Dv(ii)*zv^2*y(2,ii)+Dp(ii)*zp^2*y(3,ii)); 
PhiH=0;

res(1,ii)=-((Dh(ii+1)+Dh(ii))*(y(1,ii+1)-y(1,ii))-2*Dh(ii)*(y(1,ii)-ChH))*dxi2...
    -zh*((Dh(ii+1)*y(1,ii+1)+Dh(ii)*y(1,ii))*(y(4,ii+1)-y(4,ii))-(Dh(ii)*ChH+Dh(ii)*y(1,ii))*(y(4,ii)-PhiH))*dxi2;
res(1,ii)=zh*res(1,ii)+ze.*(-((De(ii+1)+De(ii))*(yy(1,ii+1)-yy(1,ii))-2*De(ii)*(yy(1,ii)-CeH))*dxi2...
    -ze*((De(ii+1)*yy(1,ii+1)+De(ii)*yy(1,ii))*(y(4,ii+1)-y(4,ii))-(De(ii)*CeH+De(ii)*yy(1,ii))*(y(4,ii)-PhiH))*dxi2);
%res(1,ii)=
% res(1,ii)=syp*(zh-ze*Kpeh/y(1,ii)^2).*yp(1,ii)...
%     -zh.*((Dh(ii+1)+Dh(ii)).*(y(1,ii+1)-y(1,ii))...
%     -(2*Dh(ii)).*(y(1,ii)-ChH))*dxi2...
%     -zh^2.*((Dh(ii+1).*y(1,ii+1)+Dh(ii).*y(1,ii)).*(y(4,ii+1)-y(4,ii))...
%     -(Dh(ii).*ChH+Dh(ii).*y(1,ii)).*(y(4,ii)-PhiH))*dxi2...
%     +ze*Kpeh*((De(ii+1)./(y(1,ii+1).^2)+De(ii)./(y(1,ii)).^2).*(y(1,ii+1) ...
%     -y(1,ii))-(De(ii)./(ChH).^2+De(ii)./(y(1,ii)).^2).*(y(1,ii) ...
%     -ChH))*dxi2-ze^2.*Kpeh*((De(ii+1)./y(1,ii+1) ...
%     +De(ii)./y(1,ii)).*(y(4,ii+1)-y(4,ii))...
%     -(De(ii)./ChH+De(ii)./y(1,ii)).*(y(4,ii)-PhiH))*dxi2;


res(2,ii)=syp*yp(2,ii)-((Dv(ii+1)+Dv(ii)).*(y(2,ii+1)-y(2,ii))...
    -(Dv(ii)+Dv(ii)).*(y(2,ii)-CvH))*dxi2...
    -zv*((Dv(ii+1).*y(2,ii+1)+Dv(ii).*y(2,ii)).*(y(4,ii+1)-y(4,ii))...
    -(Dv(ii).*CvH+Dv(ii).*y(2,ii)).*(y(4,ii)-PhiH))*dxi2;

res(3,ii)=syp*yp(3,ii)-((Dp(ii+1)+Dp(ii)).*(y(3,ii+1)-y(3,ii))...
    -(Dp(ii)+Dp(ii)).*(y(3,ii)-CpH))*dxi2...
    -zp*((Dp(ii+1).*y(3,ii+1)+Dp(ii).*y(3,ii)).*(y(4,ii+1)-y(4,ii))...
    -(Dp(ii).*CpH+Dp(ii).*y(3,ii)).*(y(4,ii)-CpH))*dxi2;

res(4,ii)=epsilon0*epsilonr*(y(4,ii+1)-2*y(4,ii)+PhiH).*dxi2...
     +(zh*y(1,ii)+ze*yy(1,ii)+zv*y(2,ii)+zp*y(3,ii)-0.2);

 
ii=NumGrid;

sumCF=Dh(ii)*zh*(ChOH-y(1,ii-1))+De(ii)*ze*(CeOH-yy(1,ii-1)) ...
    +Dv(ii)*zv*(CvOH-y(2,ii-1))+Dp(ii)*zp*(CpOH-y(3,ii-1));
PhiOH=y(4,ii-1)-2*dx*(ic+sumCF/2/dx)/(Dh(ii)*zh^2*y(1,ii)+De(ii)*ze^2*yy(1,ii) ...
    +Dv(ii)*zv^2*y(2,ii)+Dp(ii)*zp^2*y(3,ii)); 
%PhiOH=0.0;    
res(1,ii)=-(2*Dh(ii)*(ChOH-y(1,ii))-(Dh(ii-1)+Dh(ii))*(y(1,ii)-y(1,ii-1)))*dxi2...
    -zh*(Dh(ii)*(ChOH+y(1,ii))*(PhiOH-y(4,ii))-(Dh(ii-1)*y(1,ii-1)+Dh(ii)*y(1,ii))*(y(4,ii)-y(4,ii-1)))*dxi2;
res(1,ii)=zh*res(1,ii)+ze*(-(2*De(ii)*(CeOH-yy(1,ii))-(De(ii-1)+De(ii))*(yy(1,ii)-yy(1,ii-1)))*dxi2...
    -ze*(De(ii)*(CeOH+yy(1,ii))*(PhiOH-y(4,ii))-(De(ii-1)*yy(1,ii-1)+De(ii)*yy(1,ii))*(y(4,ii)-y(4,ii-1)))*dxi2);
% res(1,ii)=syp*(zh-ze*Kpeh/y(1,ii)^2).*yp(1,ii)...
%     -zh.*((Dh(ii)+Dh(ii)).*(CpOH-y(1,ii))...
%     -(Dh(ii-1)+Dh(ii)).*(y(1,ii)-y(1,ii-1)))*dxi2...
%     -zh^2.*((Dh(ii).*CpOH+Dh(ii).*y(1,ii)).*(PhiOH-y(4,ii))...
%     -(Dh(ii-1).*y(1,ii-1)+Dh(ii).*y(1,ii)).*(y(4,ii)-y(4,ii-1)))*dxi2...
%     +ze*Kpeh*((De(ii)./(CpOH^2)+De(ii)./(y(1,ii)^2)).*(CpOH...
%     -y(1,ii))-(De(ii-1)./(y(1,ii-1)^2)+De(ii)./(y(1,ii)^2)).*(y(1,ii)...
%     -y(1,ii-1)))*dxi2-ze^2.*Kpeh*((De(ii)./CpOH...
%     +De(ii)./y(1,ii)).*(PhiOH-y(4,ii))...
%     -(De(ii-1)./y(1,ii-1)+De(ii)./y(1,ii)).*(y(4,ii)-y(4,ii-1)))*dxi2;


res(2,ii)=syp*yp(2,ii)-((Dv(ii)+Dv(ii)).*(CvOH-y(2,ii))...
    -(Dv(ii-1)+Dv(ii)).*(y(2,ii)-y(2,ii-1)))*dxi2...
    -zv*((Dv(ii).*CvOH+Dv(ii).*y(2,ii)).*(PhiOH-y(4,ii))...
    -(Dv(ii-1).*y(2,ii-1)+Dv(ii).*y(2,ii)).*(y(4,ii)-y(4,ii-1)))*dxi2;

res(3,ii)=syp*yp(3,ii)-((Dp(ii)+Dp(ii)).*(CpOH-y(3,ii))...
    -(Dp(ii-1)+Dp(ii)).*(y(3,ii)-y(3,ii-1)))*dxi2...
    -zp*((Dp(ii).*CpOH+Dp(ii).*y(3,ii)).*(PhiOH-y(4,ii))...
    -(Dp(ii-1).*y(3,ii-1)+Dp(ii).*y(3,ii)).*(y(4,ii)-y(4,ii-1)))*dxi2;

res(4,ii)=epsilon0*epsilonr.*(PhiOH-2*y(4,ii)+y(4,ii-1)).*dxi2*2 ...
     +(zh*y(1,ii)+ze*yy(1,ii)+zv*y(2,ii)+zp*y(3,ii)-0.2);

res(1,1:NumGrid)=syp.*yp(1,1:NumGrid)+res(1,1:NumGrid).*y(1,1:NumGrid).^2./(zh.*y(1,1:NumGrid).^2-ze);
rescv(1:4*NumGrid)=0.0;

for ii=1:NumGrid
    for jj=1:4
      rescv(jj+4*(ii-1))= res(jj,ii);
  %      yp(jj,ii)=ypc(jj+5*(ii-1));
    end
end
resc=rescv(:);
end


function resI = finitial(y,icg)
global ChH CeH CpH CvH ChOH CeOH CpOH CvOH  T F R dx NumGrid ...
       Dh De Dp Dv epsilon0 epsilonr zh ze zv zp Kpeh Vm regupara
dxi2=1/(2*dx*dx); 

yy(1,1:NumGrid)=sqrt(Kpeh./(y(1,1:NumGrid).^2+regupara));
%%grid points in space
xgrid=(2:1:NumGrid-1);
resI(1,xgrid)=((Dh(xgrid+1)+Dh(xgrid)).*(y(1,xgrid+1).^2-y(1,xgrid).^2)...
    -(Dh(xgrid-1)+Dh(xgrid)).*(y(1,xgrid).^2-y(1,xgrid-1).^2))...
    +zh*((Dh(xgrid+1).*y(1,xgrid+1).^2+Dh(xgrid).*y(1,xgrid).^2).*(y(4,xgrid+1)-y(4,xgrid))...
    -(Dh(xgrid-1).*y(1,xgrid-1).^2+Dh(xgrid).*y(1,xgrid).^2).*(y(4,xgrid)-y(4,xgrid-1)));
resI(1,xgrid)=zh*resI(1,xgrid)+ze*(((De(xgrid+1)+De(xgrid)).*(yy(1,xgrid+1).^2-yy(1,xgrid).^2)...
    -(De(xgrid-1)+De(xgrid)).*(yy(1,xgrid).^2-yy(1,xgrid-1).^2))...
    +ze*((De(xgrid+1).*yy(1,xgrid+1).^2+De(xgrid).*yy(1,xgrid).^2).*(y(4,xgrid+1)-y(4,xgrid))...
    -(De(xgrid-1).*yy(1,xgrid-1).^2+De(xgrid).*yy(1,xgrid).^2).*(y(4,xgrid)-y(4,xgrid-1))));
% cfluxh(1,xgrid)=-zh.*((Dh(xgrid+1)+Dh(xgrid)).*(y(1,xgrid+1).^2-y(1,xgrid).^2)...
%     -(Dh(xgrid-1)+Dh(xgrid)).*(y(1,xgrid).^2-y(1,xgrid-1).^2))*dxi2.*y(1,xgrid).^4;
% dfluxh(1,xgrid)=-zh^2.*((Dh(xgrid+1).*y(1,xgrid+1).^2+Dh(xgrid).*y(1,xgrid).^2).*(y(4,xgrid+1)-y(4,xgrid))...
%     -(Dh(xgrid-1).*y(1,xgrid-1).^2+Dh(xgrid).*y(1,xgrid).^2).*(y(4,xgrid)-y(4,xgrid-1)))*dxi2.*y(1,xgrid).^4;
% cfluxe(1,xgrid)=ze*Kpeh*((De(xgrid+1).*(y(1,xgrid)./(y(1,xgrid+1)).^4)+De(xgrid)).*(y(1,xgrid+1).^2 ...
%     -y(1,xgrid).^2)-(De(xgrid-1).*(y(1,xgrid)./(y(1,xgrid-1)).^4)+De(xgrid)).*(y(1,xgrid).^2 ...
%     -y(1,xgrid-1).^2))*dxi2;
% dfluxe(1,xgrid)=-ze^2.*Kpeh*((De(xgrid+1).*y(1,xgrid).^2.*(y(1,xgrid)./y(1,xgrid+1)).^2 ...
%     +De(xgrid).*y(1,xgrid).^2).*(y(4,xgrid+1)-y(4,xgrid))...
%     -(De(xgrid-1).*y(1,xgrid).^2.*(y(1,xgrid)./y(1,xgrid-1)).^2+De(xgrid).*y(1,xgrid).^2).*(y(4,xgrid)-y(4,xgrid-1)))*dxi2;
% 
% resI(1,xgrid)=cfluxh(1,xgrid)+dfluxh(1,xgrid)+cfluxe(1,xgrid)+dfluxe(1,xgrid);
% resI(1,xgrid)=-zh.*((Dh(xgrid+1)+Dh(xgrid)).*(y(1,xgrid+1).^2-y(1,xgrid).^2)...
%     -(Dh(xgrid-1)+Dh(xgrid)).*(y(1,xgrid).^2-y(1,xgrid-1).^2))*dxi2...
%     -zh^2.*((Dh(xgrid+1).*y(1,xgrid+1).^2+Dh(xgrid).*y(1,xgrid).^2).*(y(4,xgrid+1)-y(4,xgrid))...
%     -(Dh(xgrid-1).*y(1,xgrid-1).^2+Dh(xgrid).*y(1,xgrid).^2).*(y(4,xgrid)-y(4,xgrid-1)))*dxi2...
%     +ze*Kpeh*((De(xgrid+1)./(y(1,xgrid+1).^4)+De(xgrid)./(y(1,xgrid).^4)).*(y(1,xgrid+1).^2 ...
%     -y(1,xgrid).^2)-(De(xgrid-1)./(y(1,xgrid-1).^4)+De(xgrid)./(y(1,xgrid).^4)).*(y(1,xgrid).^2 ...
%     -y(1,xgrid-1).^2))*dxi2-ze^2.*Kpeh*((De(xgrid+1)./y(1,xgrid+1).^2 ...
%     +De(xgrid)./y(1,xgrid).^2).*(y(4,xgrid+1)-y(4,xgrid))...
%     -(De(xgrid-1)./y(1,xgrid-1).^2+De(xgrid)./y(1,xgrid).^2).*(y(4,xgrid)-y(4,xgrid-1)))*dxi2;

%resI(1,xgrid)=resI(1,xgrid).*y(1,xgrid).^4;
% resI(2,xgrid)=((De(xgrid+1)+De(xgrid)).*(y(2,xgrid+1).^2-y(2,xgrid).^2)...
%     -(De(xgrid-1)+De(xgrid)).*(y(2,xgrid).^2-y(2,xgrid-1).^2))...
%     +ze*((De(xgrid+1).*y(2,xgrid+1).^2+De(xgrid).*y(2,xgrid).^2).*(y(5,xgrid+1).^2-y(5,xgrid).^2)...
%     -(De(xgrid-1).*y(2,xgrid-1).^2+De(xgrid).*y(2,xgrid).^2).*(y(5,xgrid).^2-y(5,xgrid-1).^2));

resI(2,xgrid)=((Dv(xgrid+1)+Dv(xgrid)).*(y(2,xgrid+1).^2-y(2,xgrid).^2)...
    -(Dv(xgrid-1)+Dv(xgrid)).*(y(2,xgrid).^2-y(2,xgrid-1).^2))...
    +zv*((Dv(xgrid+1).*y(2,xgrid+1).^2+Dv(xgrid).*y(2,xgrid).^2).*(y(4,xgrid+1)-y(4,xgrid))...
    -(Dv(xgrid-1).*y(2,xgrid-1).^2+Dv(xgrid).*y(2,xgrid).^2).*(y(4,xgrid)-y(4,xgrid-1)));

resI(3,xgrid)=((Dp(xgrid+1)+Dp(xgrid)).*(y(3,xgrid+1).^2-y(3,xgrid).^2)...
    -(Dp(xgrid-1)+Dp(xgrid)).*(y(3,xgrid).^2-y(3,xgrid-1).^2))...
    +zp*((Dp(xgrid+1).*y(3,xgrid+1).^2+Dp(xgrid).*y(3,xgrid).^2).*(y(4,xgrid+1)-y(4,xgrid))...
    -(Dp(xgrid-1).*y(3,xgrid-1).^2+Dp(xgrid).*y(3,xgrid).^2).*(y(4,xgrid)-y(4,xgrid-1)));

resI(4,xgrid)=epsilon0*epsilonr*(y(4,xgrid+1)-2*y(4,xgrid)+y(4,xgrid-1)).*dxi2*2.*y(1,xgrid).^2 ...
     +(zh*y(1,xgrid).^4+ze*Kpeh+zv*y(2,xgrid).^2.*y(1,xgrid).^2+zp*y(3,xgrid).^2.*y(1,xgrid).^2-0.2.*y(1,xgrid).^2);

%%%Set boundary condition for phi

 
%%%%boundary conditions
ii=1;
PhiH=0;
% sumCF=Dh(ii)*zh*(y(1,ii)^2-ChH^2)+De(ii)*ze*(Kpeh/y(1,ii)^2-CeH^2) ...
%      +Dv(ii)*zv*(y(2,ii)^2-CvH^2)+Dp(ii)*zp*(y(3,ii)^2-CpH^2);
% PhiH=y(4,ii)+2*dx*(icg+sumCF/2/dx)/(Dh(ii)*zh^2*y(1,ii)^2+De(ii)*ze^2*Kpeh/y(1,ii)^2 ...
%      +Dv(ii)*zv^2*y(2,ii)^2+Dp(ii)*zp^2*y(3,ii)^2); 
resI(1,ii)=-((Dh(ii+1)+Dh(ii))*(y(1,ii+1)^2-y(1,ii)^2)-2*Dh(ii)*(y(1,ii)^2-ChH^2))...
    -zh*((Dh(ii+1)*y(1,ii+1)^2+Dh(ii)*y(1,ii)^2)*(y(4,ii+1)-y(4,ii))-(Dh(ii)*ChH^2+Dh(ii)*y(1,ii)^2)*(y(4,ii)-PhiH));
resI(1,ii)=zh*resI(1,ii)+ze*(-((De(ii+1)+De(ii))*(yy(1,ii+1)^2-yy(1,ii)^2)-2*De(ii)*(yy(1,ii)^2-CeH^2))...
    -ze*((De(ii+1)*yy(1,ii+1)^2+De(ii)*yy(1,ii)^2)*(y(4,ii+1)-y(4,ii))-(De(ii)*CeH^2+De(ii)*yy(1,ii)^2)*(y(4,ii)-PhiH)));

% resI(1,ii)=-zh.*((Dh(ii+1)+Dh(ii)).*(y(1,ii+1).^2-y(1,ii).^2)...
%     -(2*Dh(ii)).*(y(1,ii).^2-ChH^2))*dxi2*y(1,ii)^4 ...
%     -zh^2.*((Dh(ii+1).*y(1,ii+1).^2+Dh(ii).*y(1,ii).^2).*(y(4,ii+1)-y(4,ii))...
%     -(Dh(ii).*ChH.^2+Dh(ii).*y(1,ii).^2).*(y(4,ii)-PhiH))*dxi2*y(1,ii)^4 ...
%     +ze*Kpeh*((De(ii+1).*y(1,ii)^4/(y(1,ii+1).^4)+De(ii)).*(y(1,ii+1).^2 ...
%     -y(1,ii).^2)-(De(ii).*y(1,ii)^4/(ChH).^4+De(ii)).*(y(1,ii).^2 ...
%     -ChH.^2))*dxi2-ze^2.*Kpeh*((De(ii+1).*y(1,ii).^4/y(1,ii+1).^2 ...
%     +De(ii).*y(1,ii).^2).*(y(4,ii+1)-y(4,ii))...
%     -(De(ii).*y(1,ii)^4/ChH.^2+De(ii).*y(1,ii).^2).*(y(4,ii)-PhiH))*dxi2;

%resI(1,ii)=resI(1,ii)*y(1,ii)^4;

resI(2,ii)=((Dv(ii+1)+Dv(ii)).*(y(2,ii+1).^2-y(2,ii).^2)...
    -(Dv(ii)+Dv(ii)).*(y(2,ii).^2-CvH.^2))...
    +zv*((Dv(ii+1).*y(2,ii+1).^2+Dv(ii).*y(2,ii).^2).*(y(4,ii+1)-y(4,ii))...
    -(Dv(ii).*CvH.^2+Dv(ii).*y(2,ii).^2).*(y(4,ii)-PhiH));

resI(3,ii)=((Dp(ii+1)+Dp(ii)).*(y(3,ii+1).^2-y(3,ii).^2)...
    -(Dp(ii)+Dp(ii)).*(y(3,ii).^2-CpH.^2))...
    +zp*((Dp(ii+1).*y(3,ii+1).^2+Dp(ii).*y(3,ii).^2).*(y(4,ii+1)-y(4,ii))...
    -(Dp(ii).*CpH.^2+Dp(ii).*y(3,ii).^2).*(y(4,ii)-PhiH));

resI(4,ii)=epsilon0*epsilonr*(y(4,ii+1)-2*y(4,ii)+PhiH).*dxi2*2.*y(1,ii)^2 ...
     +(zh*y(1,ii).^4+ze*Kpeh+zv*y(2,ii).^2*y(1,ii)^2+zp*y(3,ii).^2*y(1,ii)^2-0.2*y(1,ii)^2);

 
ii=NumGrid;

sumCF=Dh(ii)*zh*(ChOH^2-y(1,ii-1)^2)+De(ii)*ze*(CeOH^2-yy(1,ii-1)^2) ...
    +Dv(ii)*zv*(CvOH^2-y(2,ii-1)^2)+Dp(ii)*zp*(CpOH^2-y(3,ii-1)^2);
PhiOH=y(4,ii-1)-2*dx*(icg+sumCF/2/dx)/(Dh(ii)*zh^2*y(1,ii)^2+De(ii)*ze^2*yy(1,ii)^2 ...
    +Dv(ii)*zv^2*y(2,ii)^2+Dp(ii)*zp^2*y(3,ii)^2); 
%PhiOH=0.0;    
resI(1,ii)=-(2*Dh(ii)*(ChOH^2-y(1,ii)^2)-(Dh(ii-1)+Dh(ii))*(y(1,ii)^2-y(1,ii-1)^2))...
    -zh*((Dh(ii)*ChOH^2+Dh(ii)*y(1,ii)^2)*(PhiOH-y(4,ii))-(Dh(ii-1)*y(1,ii-1)^2+Dh(ii)*y(1,ii)^2)*(y(4,ii)-y(4,ii-1)));
resI(1,ii)=zh*resI(1,ii)+ze*(-(2*De(ii)*(CeOH^2-yy(1,ii)^2)-(De(ii-1)+De(ii))*(yy(1,ii)^2-yy(1,ii-1)^2))...
    -ze*((De(ii)*CeOH^2+De(ii)*yy(1,ii)^2)*(PhiOH-y(4,ii))-(De(ii-1)*yy(1,ii-1)^2+De(ii)*yy(1,ii)^2)*(y(4,ii)-y(4,ii-1))));
% resI(1,ii)=-zh.*((Dh(ii)+Dh(ii)).*(CpOH.^2-y(1,ii).^2)...
%     -(Dh(ii-1)+Dh(ii)).*(y(1,ii).^2-y(1,ii-1).^2))*dxi2*y(1,ii)^4 ...
%     -zh^2.*((Dh(ii).*CpOH.^2+Dh(ii).*y(1,ii).^2).*(PhiOH-y(4,ii))...
%     -(Dh(ii-1).*y(1,ii-1).^2+Dh(ii).*y(1,ii).^2).*(y(4,ii)-y(4,ii-1)))*dxi2*y(1,ii)^4 ...
%     +ze*Kpeh*((De(ii).*y(1,ii)^4/(CpOH.^4)+De(ii)).*(CpOH.^2 ...
%     -y(1,ii).^2)-(De(ii-1).*y(1,ii)^4/(y(1,ii-1)).^4+De(ii)).*(y(1,ii).^2 ...
%     -y(1,ii-1).^2))*dxi2-ze^2.*Kpeh*((De(ii).*y(1,ii)^4/CpOH.^2 ...
%     +De(ii).*y(1,ii).^2).*(PhiOH-y(4,ii))...
%     -(De(ii-1).*y(1,ii)^4/y(1,ii-1).^2+De(ii).*y(1,ii).^2).*(y(4,ii)-y(4,ii-1)))*dxi2;

%resI(1,ii)=resI(1,ii)*y(1,ii)^4;

resI(2,ii)=-(2*Dv(ii)*(CvOH^2-y(2,ii)^2)-(Dv(ii-1)+Dv(ii))*(y(2,ii)^2-y(2,ii-1)^2))...
    -zv*((Dv(ii)*CvOH^2+Dv(ii)*y(2,ii)^2)*(PhiOH-y(4,ii))-(Dv(ii-1)*y(2,ii-1)^2+Dv(ii)*y(2,ii)^2)*(y(4,ii)-y(4,ii-1)));

resI(3,ii)=-(2*Dp(ii)*(CpOH^2-y(3,ii)^2)-(Dp(ii-1)+Dp(ii))*(y(3,ii)^2-y(3,ii-1)^2))...
    -zp*((Dp(ii)*CpOH^2+Dp(ii)*y(3,ii)^2)*(PhiOH-y(4,ii))-(Dp(ii-1)*y(3,ii-1)^2+Dp(ii)*y(3,ii)^2)*(y(4,ii)-y(4,ii-1)));

resI(4,ii)=epsilon0*epsilonr*(y(4,ii-1)-2*y(4,ii)+PhiOH).*dxi2*2.*y(1,ii)^2 ...
     +(zh*y(1,ii).^4+ze*Kpeh+zv*y(2,ii).^2*y(1,ii)^2+zp*y(3,ii).^2*y(1,ii)^2-0.2*y(1,ii)^2);
end
