function molality=methane_solubility(P,T,NaCl)

c1= 4.30310345E+1;
c2=-6.83277221E-2;
c3=-5.68718730E+3;
c4= 3.56636821E-5;
c5=-5.79133791E+1;
c6= 6.11616662E-3;
c7=-7.85528103E-4;
c8=-9.42540759E-2;
c9= 1.92132040E-2;
c10= -9.17186899E-06;

lambda_c1=9.92230792E-2;
lambda_c2=2.57906811E-5;
lambda_c3=0.0;
lambda_c4=0.0;
lambda_c5=0.0;
lambda_c6=0.0;
lambda_c7=0.0;
lambda_c8=1.83451402E-2;
lambda_c9=0.0;
lambda_c10=-8.07196716E-6;

xi_c1=-6.23943799E-3;
xi_c2=0.0;
xi_c3=0.0;
xi_c4=0.0;
xi_c5=0.0;
xi_c6=0.0;
xi_c7=0.0;
xi_c8=0.0;
xi_c9=0.0;
xi_c10=0.0;

P=P*10;
T=273.15+T;

% P is in bar and T is in K. 

Tc=190.6;
Pc=46.41;
Pr=P/Pc;
Tr=T/Tc;

par=c1+c2*T+c3/T+c4*T^2+c5/(680-T)+...
    c6*P+c7*P*log(T)+c8*P/T+c9*P/(680-T)+c10*P^2/T;

lambda_par=lambda_c1+lambda_c2*T+lambda_c3/T+lambda_c4*T^2+lambda_c5/(680-T)+...
    lambda_c6*P+lambda_c7*P*log(T)+lambda_c8*P/T+lambda_c9*P/(680-T)+lambda_c10*P^2/T;

xi_par=xi_c1+xi_c2*T+xi_c3/T+xi_c4*T^2+xi_c5/(680-T)+...
    xi_c6*P+xi_c7*P*log(T)+xi_c8*P/T+xi_c9*P/(680-T)+xi_c10*P^2/T;

%Vr=fsolve('duan',Tr/Pr,optimset('Display','off'),Pr,Tr);
Vr=calc_Vr(Pr,Tr);
Z=Pr*Vr/Tr;


a1= 8.72553928E-2;
a2=-7.52599476E-1;
a3= 3.75419887E-1;
a4= 1.07291342E-2;
a5= 5.49626360E-3;
a6=-1.84772802E-2;
a7= 3.18993183E-4;
a8= 2.11079375E-4;
a9= 2.01682801E-5;
a10=-1.65606189E-5;
a11= 1.19614546E-4;
a12=-1.08087289E-4;
alpha= 4.48262295E-2;
beta= 7.53970000E-1;
gamma= 7.7167000E-2;

B=a1+a2/Tr^2+a3/Tr^3;
C=a4+a5/Tr^2+a6/Tr^3;
D=a7+a8/Tr^2+a9/Tr^3;
E=a10+a11/Tr^2+a12/Tr^3;
F=alpha/Tr^3;
G=F/(2*gamma)*(beta+1-(beta+1+gamma/Vr^2)*exp(-gamma/Vr^2));

fugacity_coef=exp(Z-1-log(Z)+B/Vr+C/(2*Vr^2)+D/(4*Vr^4)+E/(5*Vr^5)+G);
%fugacity=fugacity_coef*P;

molality=P/exp(par+2*lambda_par*NaCl+xi_par*NaCl^2-log(fugacity_coef));
    



