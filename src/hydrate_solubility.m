function molality=hydrate_solubility(P0,T0,NaCl) 
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

P=P0*10;
T=273.15+T0;

lambda_par=lambda_c1+lambda_c2*T+lambda_c3/T+lambda_c4*T^2+lambda_c5/(680-T)+...
    lambda_c6*P+lambda_c7*P*log(T)+lambda_c8*P/T+lambda_c9*P/(680-T)+lambda_c10*P^2/T;

xi_par=xi_c1+xi_c2*T+xi_c3/T+xi_c4*T^2+xi_c5/(680-T)+...
    xi_c6*P+xi_c7*P*log(T)+xi_c8*P/T+xi_c9*P/(680-T)+xi_c10*P^2/T;

% P is in bar and T is in K. 

par=c1+c2*T+c3/T+c4*T^2+c5/(680-T)+...
    c6*P+c7*P*log(T)+c8*P/T+c9*P/(680-T)+c10*P^2/T;

aw=water_activity(NaCl);
fugacity=calc_fugacity(P0,T0,aw);
molality=fugacity*exp(-par-2*lambda_par*NaCl-xi_par*NaCl^2);





