function y=fugacity_methane(x,P0,T0,aw)

R=83.1441;

% vi is the number of cavitis of type i per molecular of water
% 46 water, 2 small cavity and 6 large cavity per unit cell.
v1=2/46;
v2=6/46;

% Small cavity 
A1=0.7228e-3;   % K/atm
B1=3187;   % K

% large cavity
A2=23.35e-3;  % K/atm
B2=2653;   % K


delta_mu0=12640;
delta_H0=-48580;
delta_V0=4.6;
delta_Cp=-391.6;

P=P0*10;
T=273.15+T0;
T_bar=(T+273.15)/2;


left=delta_mu0/(R*273.15) - ...
     (delta_H0-delta_Cp*273.15)/R*(-1/T+1/273.15) - ...
     delta_Cp/R*(log(T)-log(273.15)) + ...
     delta_V0/(R*T_bar)*P - log(aw); % + ...

 C1=A1/T*exp(B1/T);
 C2=A2/T*exp(B2/T);
 
 y=left+v1*log(1-C1*x/(1+C1*x))+v2*log(1-C2*x/(1+C2*x));
 
 









