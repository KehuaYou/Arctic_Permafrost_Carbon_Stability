function rho=brine_density(P,T,m,xx)
phi_V0=18.22;
A=2.5830;
B=-0.2624;
C=0.01598;
D=1.4270;
E=-0.4265;
F=-0.004437;
G=0.001353;
H0=-631.89E-4;
H1=3.254E-4;
H2=0.0047E-4;
I0=470.26E-4;
I1=-2.668E-4;
I2=0.0039E-4;

Av=1.8743;

term1=0;
term2=0;
T0=298;
P0=1.01325;

m=m/2;
P=P*10;
T=273.15+T;

%m=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Calc rhow for fresh water                                         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_P=P/10*1e6/(9.81E4);
AT= 5.916365 - 0.01035794*T + 0.9270048E-5*T^2 ...
   -1127.522/T + 100674.1/T^2;
BT= 0.5204914E-2 - 0.10482101E-4*T + 0.8328532E-8*T^2 ...
   -1.1702939/T + 102.2783/T^2;
CT= 0.118547E-7-0.6599143E-10*T;
v0=AT-new_P*BT-new_P^2*CT;
rhow=1/v0; % gm/cm^3

%rhow=1;

T_interp=[T0:-1:T];
dT=1;
[row,n]=size(T_interp);
for i=1:(n-1)
    T_temp=(T_interp(i)+T_interp(i+1))/2;
    term1=term1+(D+E*sqrt(m))+(F+G*sqrt(m))*T_temp;
end
H=H0+H1*T+H2*T^2;
I=I0+I1*T+I2*T^2;
term2=(H+I*sqrt(m))*(P-1);

phi_Vm=phi_V0+3*sqrt(3)*Av*sqrt(m)/(1+sqrt(3)*sqrt(m))+A*m+B*m^2+C*m^3;
phi_V=phi_Vm+term1-term2;

Mc=110.99; % formula weight of CaCl2 (gm/mol)
Mch4=16;
phi_Vch4=3.7E-05*1e6;
%phi_Vch4=0;
%mch4=0;
%m=0;
rho_temp=(m*Mc+1000)/(m*phi_V+1000/rhow);
m_total=1000+m*Mc+xx*Mch4;
V_total=(1000+m*Mc)/rho_temp+xx*phi_Vch4;
%rho=(m*Mc+1000+xx*Mch4)/(m*phi_V+1000/rhow+xx*phi_Vch4)*1000;
rho=m_total/V_total*1000;


