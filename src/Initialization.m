
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global COL ROW
global INDC1 INDC2 INDC3
global p_salinity T_salinity salinity
global dpth dx dy dz  vb
global cnstx cnstTx
global pw pw_0 
global T T_0
global sw sw_0 sh sh_0 sg sg_0 si si_0
global cl cl_0 cm cm_0
global dnw dnw_0 dng dng_0 dnh dns dni g
global Mh2o Mch4
global phi phi_i  phi_0
global kx ky k_temp
global krw krg sgr swr wn gn 
global Dmethanex Dsaltx
global vsw vsw_0 vsg vsg_0 viscosity_g viscosity_w
global pcgw pcgw_0 Pd0 
global tsg tpcgw nrpermg
global lambda lambda_g lambda_h lambda_w lambda_s lambda_i
global Cp_g Cp_h Cp_w Cp_s Cp_i
global qg_biogenic C_stable C_labile
global qw qg qt
global pw_rate T_top cm_top cl_top
global eps
global pw_ini T_ini pw_scale T_scale cl_scale cm_scale



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Phase Diagram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_salinity_permafrost.mat   % calculated using the program "Kehua_phase_curve.m"
T_salinity=(-20:0.2:20);
salinity=(0:0.5:22)/100;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------- density & mass fraction  ------------
g=9.81;     % gravitational constant
dnh=912;    % density of hydrate
dni=917;    % ice density at 0 degree c
dns=2750;   % density of solid grain
Mh2o=5.75*0.018/(0.016+5.75*0.018); % molar weight fraction of water in hydrate
Mch4=0.016/(0.016+5.75*0.018);      % molar weight fraction of methane in hydrate

% ---------- Flow properties ------------
sgr=0.02;   % residual gas saturation
swr=0.1;    % residual water saturation
wn=4;       % constant used in equations to calculate relative water permeability 4 originally
gn=2;       % constant used in equations to calculate relative gas permeability
viscosity_g = 2.0e-5;   % dynamic viscosity of gas
viscosity_w = 1.31e-3;  % dynamic viscosity of water

% ---------- Gas-Liquid Water capillary pressure------------
nrpermg=12;
tsg=[0 0.02 0.04 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
tpcgw=[0.20 0.21 0.22 0.23 0.25 0.3 0.4 0.5 0.65 0.8 0.95 1.15]*5;

%------------ thermal constant ----------------------
Cp_g=3500; % J/K*kg - heat capacity of methane gas
Cp_h=2100; % J/K*kg - heat capacity of hydrate
Cp_w=4200; % J/K*kg - heat capacity of water
Cp_s=730;%1381; % J/K*kg - heat capacity of quartz
Cp_i=2108;

lambda_g=0.03281; % W/m*K - thermal conductivity of gas
lambda_h=0.49; % W/m*K - thermal conductivity of hydrate
lambda_w=0.58; % W/m*K - thermal conductivity of water
lambda_s=2.3;%3.1;%1.78;%3.0;%5.5;%3.0;%1.60; % W/m*K - thermal conductivity of quartz
lambda_i=2.2;

%------------ numerical method  ----------------------
eps=1E-6;   % change in parameter when calculate derivative and the Jacobian matrix

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qw=0;           % Bottom water flux in kg/m2/yr
qg=0;           % Bottom gas flux in kg/m2/yr 2e-6 m/s
qt=60.4e-3;     % Geothermal heat flux in W/m2


T_top=-1.3;     % Seafloor temperature 
pw_rate=120/(18e3*86400*365);   % Seafloor hydrostatic pore pressure increasing rate
cm_top=0;       % Seafloor dissolved methane concentration
cl_top=0.032;   % Seafloor salinity

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Geometery
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ROW=1;      % number of horizontal grids
COL=121;    % number of vertical grids
dx=10.*ones(COL,ROW);   % vertical grid size
dy=10*ones(COL,ROW);    % horizontal grid size
dz=1*ones(COL,ROW);     % horizontal grid size

seafloor=0;         % water depth at initial condition
model_top=0;        % depth of mdoel top

dpth_temp=zeros(COL,ROW);
dpth_temp(1,:)=model_top;
for i=1:COL-1
    dpth_temp(i+1,:)=dpth_temp(i,:)+(dx(i,:)+dx(i+1,:))*0.5;
end
dpth=dpth_temp;     % depth of each grid


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Initialization 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------initialize saturation--------------------------------
for i=1:COL
    for j=1:ROW
        sw_0(i,j)=1;
        sh_0(i,j)=0;
        sg_0(i,j)=0;
        si_0(i,j)=0;
    end
end

%-----------------------initial temeprature--------------------------------
T_surf=-20;         % Ground surface temperature in degC
cl_bp = 0.032;      % Salinity at the base of permafrost
dpth_bp = 650;      % Thickness of permafrost
col_bp = dpth_bp/dx(1,1)+1; % Grid number at the base of permafrost
Tfz_bp=-cl_bp*(164.40*cl_bp+49.462);    % Temperature at the base of permafrost
T_grad1=(Tfz_bp-T_surf)/dpth_bp;        % Temperature gradient within permafrost in degC/m
T_grad2=0.0356;                         % Temperature gradient below permafrost

for i=1:COL
    if dpth(i,1)<=dpth_bp
        T(i,1)=T_surf+T_grad1*dpth(i,1);
    else
        T(i,1)=T(col_bp,1)+T_grad2*(dpth(i,1)-dpth_bp);
    end
end
T_0=T;

%--------------------------initial salinity-------------------------------
cl_0=0.032.*ones(COL,ROW);      % Initial salinity below permafrost

for i=1:col_bp-1 
    si_0(i,1) = 0.8;        % Initial ice saturation
    cl_temp=((49.462^2-4*164.49*T_0(i,1))^0.5-49.462)/(2*164.49);
    cl_0(i,1)=cl_temp;      % Initial salinity within permafrost
end


sh_0(21:23,1)=0.5;      % Initial methane hydrate saturation
sh_0(66:68,1)=0.5;
sh_0(103:105,1)=0.5;

si_0(21:23,1) = si_0(21:23,1) - sh_0(21:23,1); % Modified initial ice saturation
sw_0=1-si_0-sh_0-sg_0;                          % Initial water saturation

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Sediment hydraulic properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sediment_perm=1e-15;                                        % permeability in m-2
phi_i=0.35;                                                 % Sediment porosity
Pd0=72e-3*0.15*(phi_i/sediment_perm)^0.5*ones(COL,ROW);     % Capillary-entry pressure in Pascal, the J-function
k_temp=sediment_perm*ones(COL,ROW);                         % Intrinsic permeability 

%----------------------initial pressure & capillary pressure---------------------------
for i=1:COL
    for j=1:ROW
        pw_0(i,j)=0.1e6+dpth(i,j)*1030*g; %initial hydrostatic pore pressure
        
        if sh_0(i,j)>0 || si_0(i,j)>0
            tmp_ratio=sqrt(1/(1+si_0(i,j)+sh_0(i,j)+2*(1-si_0(i,j)-si_0(i,j))/log(si_0(i,j)+sh_0(i,j))));
        else
            tmp_ratio=1;
        end
        if tmp_ratio==0
            pcgw_0(i,j)=0;
        else
            pcgw_0(i,j)=Pd0(i,j)*tmp_ratio*interp(tsg,tpcgw,nrpermg,sg_0(i,j)/(1-sh_0(i,j)-si_0(i,j)),1); % initial capillary pressure
        end
    end
end

%---------------initial dissolved methane concentration, density, viscosity and porosity------------------------------
for i=1:COL
    for j=1:ROW
        if sh_0(i,j)+si_0(i,j)==1
            cm_0(i,j)=0;                % Initial dissolved methane concentration
        else
            if sh_0(i,j)>0
                cm_0(i,j)=hydrate_solubility(pw_0(i,j)/1e6,T_0(i,j),2*cl_0(i,j)/0.11099/(1-cl_0(i,j)))*0.016;
            elseif sg_0(i,j)>0
                cm_0(i,j)=methane_solubility(pw_0(i,j)/1e6,T_0(i,j),2*cl_0(i,j)/0.11099/(1-cl_0(i,j)))*0.016;
            else
                cm_0(i,j)=0;
            end
        end
        
        dnw_0(i,j)=brine_density(pw_0(i,j)/1e6,T_0(i,j),2*cl_0(i,j)/0.11099/(1-cl_0(i,j)), cm_0(i,j)/0.016);     % initial brine density
        vsw_0(i,j)=viscosity_w;                             % water viscosity
        dng_0(i,j)=gas_density(pw_0(i,j)/1e6,T_0(i,j));     % initial gas density
        vsg_0(i,j)=viscosity_g;                             % gas viscosity
        phi_0(i,j)=phi_i;                                   % porosity
    end
end
%-----------------initial sediment effective permeability------------------------
for i=1:COL
    for j=1:ROW
        vb(i,j)=dx(i,j)*dy(i,j)*dz(i,j);    % Volume of eacg grid block
        if sh_0(i,j)>0 || si_0(i,j)>0
            tmp=(1-si_0(i,j)-sh_0(i,j))^2;
        else
            tmp=1;
        end
        kx(i,j)=k_temp(i,j)*tmp;
        ky(i,j)=k_temp(i,j)*tmp;
    end
end

for i=1:COL
    for j=1:ROW
        if i < COL
            cnstx(i,j)=2.0 * dz(i+1,j)*dy(i+1,j)*dz(i,j)*dy(i,j)*kx(i+1,j) * kx(i,j)/ ...
                (dx(i,j)*dz(i+1,j)*dy(i+1,j)*kx(i+1,j) + dx(i+1,j)*dz(i,j)*dy(i,j)*kx(i,j));
            cnstTx(i,j)=2.0*dz(i+1,j)*dy(i+1,j)*dz(i,j)*dy(i,j)/ ...
                (dx(i,j)*dz(i+1,j)*dy(i+1,j) + dx(i+1,j)*dz(i,j)*dy(i,j));
        else
            cnstx(i,j)=0.0;
            cnstTx(i,j)=0.0;
        end
    end
end
%---------------------- diffusion coefficient ----------------------------
D0=3.6e-10.*ones(COL,ROW);      % Effective diffusion coefficient in pore water
Dmethanex=zeros(COL,ROW);
Dsaltx=zeros(COL,ROW);
for i=1:COL-1
        Dmethanex(i,1)=dy(i,1)*dz(i,1)/((dx(i,1)+dx(i+1,1))*0.5)*D0(i,1);
        Dsaltx(i,j)=dy(i,1)*dz(i,1)/((dx(i,1)+dx(i+1,1))*0.5)*D0(i,1);
end
Dmethanex(COL,1)=Dmethanex(COL-1,1);
Dsaltx(COL,1)=Dsaltx(COL-1,1);

%-----------------initial methanogeneis & organic carbon content-----------
qg_biogenic=zeros(COL,1);               % Initial methane production rate
C_total = 0.035.*(dpth<=50) + 0.01.*(dpth>50);  % Initial total organic carbon content
f_labile = 0.005;                       % Fraction of labile component 
C_labile = C_total*f_labile;            % Initial labile organic carbon content
C_stable = C_total*(1-f_labile);        % Initial stable organic carbon content

%------------------------------------------------------------------------
pw=pw_0;
cm=cm_0;
cl=cl_0;
dnw=dnw_0;
vsw=vsw_0;
dng=dng_0;
vsg=vsg_0;
phi=phi_0;
sw=sw_0;
sh=sh_0;
sg=sg_0;
si=si_0;
pcgw=pcgw_0;

%-----------------initial water & gas relative permeability----------------
for i=1:COL
    for j=1:ROW
        if sh(i,j)+si(i,j)==1
            krg(i,j)=0;
            krw(i,j)=0;
        else
            if sg(i,j)/(1-sh(i,j)-si(i,j))>sgr
                krg(i,j)=(sg(i,j)/(1-sh(i,j)-si(i,j))-sgr)^gn;     % initial gas phase relative permeability
            else
                krg(i,j)=0;
            end
            
            if sw(i,j)/(1-sh(i,j)-si(i,j))>swr
                krw(i,j)=(sw(i,j)/(1-sh(i,j)-si(i,j))-swr)^wn;      % initial liquid water relative permeability
            else
                krw(i,j)=0;
            end
        end
        
    end
end

%----------------initial bulk thermal conductivity, and phase mobility -----------------------------
lambda=(1-phi)*lambda_s+phi.*(sw*lambda_w+sg*lambda_g+sh*lambda_h+si*lambda_i); %overall thermal conductivity of the porous media
for i=1:COL
    for j=1:ROW
        calc_trans(i,j);    % Calc_trans = calculate phase mobility: gas, liquid water
    end
end

% ---------initial index (INDC1, INDC2, INDC3) values----------------------
for i=1:COL
    for j=1:ROW
        % INDC3=1, water stable zone; INDC3=2, ice stable zone; INDC3=4,
        % water and ice phase boundary; INDC3=3, transfer from water/ice
        % stable zone to ice water phase boundary
        
        if sw(i,j)>0 && si(i,j)==0
            INDC3(i,j)=1;
        elseif sw(i,j)>0 && si(i,j)>0
            INDC3(i,j)=4;
        elseif sw(i,j)==0 && si(i,j)>0
            INDC3(i,j)=2;
        end
        
        % INDC1=1, hydrate stable zone; INDC1=0, gas stable zone; phase
        % boundary, INDC1=0.
        if INDC3(i,j)==2
            if pw(i,j)/1e6 > interpolation_pressure(p_salinity,T_salinity,salinity,cl(i,j),T(i,j))
                INDC1(i,j)=1;
            else
                INDC1(i,j)=0;
            end
            
        else
            if cl(i,j) < interpolation2(p_salinity,T_salinity,salinity,pw(i,j)/1e6,T(i,j))
                INDC1(i,j)=1;
            else
                INDC1(i,j)=0;
            end
        end
        
        % INDC2=1, one phase present; INDC2=3, two phases present; INDC2=5,
        % three phases present; INDC2=7, four phases present.
        INDC2(i,j)=1;
        if (sw(i,j)>0)+(sg(i,j)>0)+(sh(i,j)>0)+(si(i,j)>0) ==2
            INDC2(i,j)=3;
        end
        if (sw(i,j)>0)+(sg(i,j)>0)+(sh(i,j)>0)+(si(i,j)>0) ==3
            INDC2(i,j)=5;
        end
        if (sw(i,j)>0)+(sg(i,j)>0)+(sh(i,j)>0)+(si(i,j)>0) ==4
            INDC2(i,j)=7;
        end
    end
end

%-----------initial total H2O, CH4 and salt masses-------------------------
init_w=sum(vb.*phi.*dnw.*sw.*(1-cl)+vb.*phi.*dni.*si+vb.*phi.*dnh.*sh*Mh2o);     % initial total water mass M/L^3
init_g=sum(vb.*phi.*dnh.*sh.*Mch4) + sum(vb.*phi.*dnw.*sw.*(1-cl).*cm) + sum(vb.*phi.*dng.*sg);      % initial total gas mass term
init_c=sum(vb.*phi.*dnw.*sw.*cl);       % initial total hydrate mass

%------------------- For non-dimensionlization -------------------------
pw_ini=pw;
T_ini=T;
pw_scale=1010*g*500;
T_scale=20;
cl_scale=0.03;
cm_scale=0.001;

% ------------- Top Boundary --------------------------
T(1)=T_top;
cm(1)=cm_top;
cl(1)=cl_top;
sw(1) = 1;
si(1) = 0;

% ------------- Display Initial Condition --------------------------
Initial_Results = ['depth   '' pw ' '  T '   ' Sw '  '   si' '   Sg'  '   sh' '   cl'   '     cm'  '   INDC1'    '    INDC2'  '   INDC3']
Initial_Results = [dpth(:,1)./1e3  pw(:,1)./1e6  T(:,1)   sw(:,1)  si(:,1)   sh(:,1) cl(:,1)   cm(:,1) C_labile(:,1)  C_stable(:, 1) INDC1(:,1)   INDC2(:,1)   INDC3(:,1)]
