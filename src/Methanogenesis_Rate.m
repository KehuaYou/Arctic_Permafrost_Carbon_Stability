global COL T_old
global C_labile C_stable C_labile_old C_stable_old
global qg_biogenic
global dt dns phi
global dpth

%% Introductory Carbon Balance Model
C_input = 0; % organic carbon input
k1= 0.3e-7; % per second
k2 = k1/0.5e4;% k1*1e-2; % per second
h=0.7;  % 40% of labile organic carbon humidification into stable pool

%% Temperature response factor
Tmin = -17;
Tref = 4;
nt = 4;
%%
for i = 1:COL
    if dpth(i,1) <= 50
        r_temperature = ((T_old(i)-Tmin)./(Tref-Tmin))^nt;
    else
        r_temperature = 1e-4*((T_old(i)-Tmin)./(Tref-Tmin))^nt;
    end

    qg_biogenic(i,1) = 0.5*16/12*dns*(1-phi(i,1))*((1-h)*k1*r_temperature*C_labile_old(i,1) + k2*r_temperature*C_stable_old(i,1));
    C_labile(i,1) = C_labile_old(i,1)*exp(-k1*r_temperature*dt);
    C_stable(i,1) = (C_stable_old(i,1)-h*k1*r_temperature*C_labile_old(i,1)/(r_temperature*(k2-k1)))*exp(-k2*r_temperature*dt) + ...
        h*k1*r_temperature*C_labile_old(i,1)/(r_temperature*(k2-k1))*exp(-k1*r_temperature*dt);
end
