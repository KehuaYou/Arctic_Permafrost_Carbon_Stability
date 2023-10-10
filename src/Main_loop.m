clc
clear all
global COL ROW
global INDC1 INDC2 INDC3 INDC2_old
global p_salinity T_salinity salinity
global dpth dy dz  
global pw pw_0
global T T_0
global sw sw_0 sh sh_0 sg sg_0 si si_0
global cl cl_0 cm cm_0
global dnw dnw_0 dng dng_0 g
global phi phi_0
global kx ky
global krw krg
global Dmethanex
global vsw vsw_0 vsg vsg_0
global pcgw pcgw_0
global lambda
global C_stable C_labile
global pw_rate T_top cm_top
global dt t_flag cycle go_back
global pw_top
global tgx tgx1 twx  twx1

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%      Load  Initializatoin File             %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option-1: run from initialization file
Initialization;                     % load initial file
string = ('t0kyr.mat');             % save initial input
save(string)
timestep=1;                         % the start timestep


% Option-2: run from saved .mat file
% load('t19.2kyr.mat')              % load initial input from .mat file
% timestep=timestep+1;              % the start timestep=current timestep+1
% string_transient=('19.2kyr.mat');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%      time discretizaton            %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt_ini=86400*365*0.5;             % Initial time step
dt_cut_inner=0.1;               % when cut time step, dt=dt*dt_cut
dt_cut_outter=0.1;              % when we go back to previous saved timestep, dt_ini=dt_ini*dt_cut_outter
dt_minimum=1;                   % the smallest timestep is 0.1 year
N=100000;                       % total number of timestep we will run
N_save=20;                      % we save our results every N_save timestep

%%%%% you need to comment this two lines if load from nonzero time
time=zeros(N+1,1);          % initializing time for each step
time(1)=0;                  % we start at time=zero

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   Main Subroutine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_iter=3;         % maximum number of iterations for new-raphson to converge
go_back=0;          % index of primariy variable shift or not
break_flag=0;

% ------------- Time-dependent Results to save -----------------
flux_advection = zeros(N+1,1);      % seafloor methane emission rate by advection with water flow
flux_diffusion=zeros(N+1,1);        % seafloor methane emission rate by diffusion in water
flux_flow=zeros(N+1,1);             % seafloor methane emission rate by free gas flow
flux_total=zeros(N+1,1);            % total methane flux at the seafloor
accumulative_advection=zeros(N+1,1);    % cumulative methane emission by advection with water flow
accumulative_diffusion=zeros(N+1,1);    % cumulative methane emission by diffusion in water
accumulative_flow=zeros(N+1,1);         % cumulative methane emission by free gas flow
accumulative_total=zeros(N+1,1);        % total cumulative methane emission at the seafloor


while timestep<=N

    disp(timestep)
    iteration=0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%       Save the data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pw_old=pw;
    cm_old=cm;
    cl_old=cl;
    sw_old=sw;
    sh_old=sh;
    sg_old=sg;
    si_old=si;
    pcgw_old=pcgw;
    dnw_old=dnw;
    dng_old=dng;
    krw_old=krw;
    krg_old=krg;
    kx_old=kx;
    ky_old=ky;
    T_old=T;
    C_labile_old=C_labile;
    C_stable_old=C_stable;
    INDC1_old=INDC1;
    INDC2_old=INDC2;
    INDC3_old=INDC3;
    lambda_old=lambda;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % at new time step, we return to original dt
    dt=dt_ini;
    time(timestep+1)=time(timestep)+dt;

    cycle=0;
    loop=1;
    while loop==1

        go_back=0;
        iteration=iteration+1;

        %*************** Top Boundary Condition ****************************
        pw_top=0.1e6+pw_rate*time(timestep)*1010*g;
        pw(1)=pw_top;
        %-------- Future seabed temperature, SSP2-4.5 here----------------
        %----- You need to edit here to reflect other warming senario-----
        if time(timestep) >= 86400*365*18e3 && time(timestep)<86400*365*18.1e3
            T(1) = T_top + 2/(86400*365*100) * (time(timestep)-86400*365*18e3);
        elseif time(timestep)>=86400*365*18.1e3
            T(1) = T_top + 2;
        end

        %************** Microbial methane generation ***********************
        Methanogenesis_Rate;

        %************** Solve the PDEs *************************************
        newton_raphson;

        %************** Check if Jacobian matrix blows up **************
        nan_flag=0;
        for i=1:COL
            if isnan(pw(i,1)) && INDC2(i,1)<6  % if the Jacobian matrix blows up, giving NAN values
                nan_flag=1;
                go_back=go_back+1;
                iteration=max_iter+1;   % get out of the loop
                t_flag=1;               % decrease the timestep
            end

            if isnan(T(i,1)) && INDC2(i,1)>=6  % if the Jacobian matrix blows up, giving NAN values
                nan_flag=1;
                go_back=go_back+1;
                iteration=max_iter+1;   % get out of the loop
                t_flag=1;               % decrease the timestep
            end
        end

        %************** Update Parameters **********************************
        if nan_flag==0
            ppt_update;
        end

        %************** Primary Variable Shift *******************
        if iteration > max_iter
            if nan_flag==0
                cycle=cycle+1;
                Phase_Zone_Judge;
                Management;
            end

            %************** Re-iteration or Go to next Time step **********
            if dt>=dt_minimum
                if t_flag==1
                    dt=dt_cut_inner*dt;
                    display(dt)
                    time(timestep+1)=time(timestep)+dt;
                    INDC1=INDC1_old;
                    INDC2=INDC2_old;
                    INDC3=INDC3_old;
                    t_flag=0;
                    go_back=go_back+1;
                    cycle=0;
                end

                if go_back>0
                    pw=pw_old;
                    cm=cm_old;
                    cl=cl_old;
                    sw=sw_old;
                    sh=sh_old;
                    sg=sg_old;
                    si=si_old;
                    krw=krw_old;
                    krg=krg_old;
                    kx=kx_old;
                    ky=ky_old;
                    pcgw=pcgw_old;
                    dnw=dnw_old;
                    dng=dng_old;
                    T=T_old;
                    lambda=lambda_old;
                    C_labile=C_labile_old;
                    C_stable=C_stable_old;

                    for i=1:COL
                        for j=1:ROW
                            calc_trans(i,j);
                        end
                    end
                    pw_0=pw;
                    cm_0=cm;
                    cl_0=cl;
                    sw_0=sw;
                    sh_0=sh;
                    sg_0=sg;
                    si_0=si;
                    pcgw_0=pcgw;
                    phi_0=phi;
                    dnw_0=dnw;
                    vsw_0=vsw;
                    dng_0=dng;
                    vsg_0=vsg;
                    T_0=T;
                    loop=1;
                    iteration=0;
                else
                    loop=0;
                    pw_0=pw;
                    cm_0=cm;
                    cl_0=cl;
                    sw_0=sw;
                    sh_0=sh;
                    sg_0=sg;
                    si_0=si;
                    pcgw_0=pcgw;
                    dnw_0=dnw;
                    dng_0=dng;
                    T_0=T;

                    for i=1:COL
                        for j=1:ROW
                            calc_trans(i,j);
                        end
                    end
                end % end of if go_back==1

            else

                break

            end

        end
    end

    if loop==0
        %****************** Decide Stable Phase zone **********************
        for i=1:COL
            for j=1:ROW
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
            end
        end

        %************** Calcualted Seabed Methane Emission ****************
        tmp_w=twx(1,1)*(pw(2,1)-pw_top) - twx1(1,1)*(dpth(2,1)-dpth(1,1));
        tmp_g=tgx(1,1)*(pw(2,1)+pcgw(2,1)-pw_top) - tgx1(1,1)*(dpth(2,1)-dpth(1,1));

        flux_advection(timestep+1) = (tmp_w>0)*tmp_w*dnw(2,1)*(1-cl(2,1))*cm(2,1)/(dy(1,1)*dz(1,1));
        flux_diffusion(timestep+1) = Dmethanex(1,1)*(phi(1,1)*sw(1,1)+phi(2,1)*sw(2,1))/2*(dnw(1,1)+dnw(2,1))/2*(cm(2,1)-cm_top)/(dy(1,1)*dz(1,1));
        flux_flow(timestep+1) = (tmp_g>0)*tmp_g*dng(2,1)/(dy(1,1)*dz(1,1));
        flux_total(timestep+1)=flux_advection(timestep+1) + flux_diffusion(timestep+1) + flux_flow(timestep+1);


        accumulative_advection(timestep+1)=accumulative_advection(timestep)+flux_advection(timestep+1)*dt;
        accumulative_diffusion(timestep+1)=accumulative_diffusion(timestep)+flux_diffusion(timestep+1)*dt;
        accumulative_flow(timestep+1)=accumulative_flow(timestep)+flux_flow(timestep+1)*dt;
        accumulative_total(timestep+1)=accumulative_total(timestep)+flux_total(timestep+1)*dt;
        
        %*********************** Save & Display Results *******************
        if ~mod(timestep,N_save)
            if break_flag==1
                break_flag=0;
                dt_ini=dt_ini/dt_cut_outter;
            end
            string = ['t',num2str(time(timestep+1)/(86400*365*1e3)),'kyr','.mat'];
            string_transient=string;
            save(string)

            % --------------- Display Results --------------------------
            junk=['depth   '' pw ' '  T '   ' Sw '  '   si' '   Sg'  '   sh' '   cl'   '     cm'  '   INDC1'    '    INDC2'  '   INDC3']
            junk=[dpth(:,1)./1e3 pw(:,1)./1e6  T(:,1)   sw(:,1)  si(:,1) sg(:,1)  sh(:,1) cl(:,1)   cm(:,1)   INDC1(:,1)   INDC2(:,1)   INDC3(:,1) C_stable(:,1)]

        end
        timestep=timestep+1;
    end

    %*********************** Go Back to Previously Saved Step *************
    if dt<dt_minimum
        if timestep<=N_save
            load('initial_input.mat');
            timestep=1;
        else
            load(string_transient);
        end
        disp('we have gone back to previous saved timestep')
        dt_ini=dt_ini*dt_cut_outter;
        break_flag=1;

        if timestep>1
            timestep=timestep+1;
        end
    end

end














