function aw=water_activity(molality)

A_phi=0.3915; % at temperature of 25 deg C
b=1.2; % unit is kg1/2 / mol1/2

m=molality;

mu_M=1;
mu_X=1;
mu=2;

z_M=+1;
z_X=-1;

I=1/2*(m*z_M^2+m*z_X^2);

f_phi=-A_phi*sqrt(I)/(1+b*sqrt(I));

alpha=2.0; % unit is kg1/2 / mol1/2
B0=0.0765;
B1=0.2664;

B_MX=B0+B1*exp(-alpha*sqrt(I));

C_MX=0.00127;

phi=1+abs(z_M*z_X)*f_phi+m*(2*mu_M*mu_X/mu)*B_MX+m^2*(2*(sqrt(mu_M*mu_X))^3/mu)*C_MX;

aw=exp(-18*mu*m/1000*phi);