%function rrr=newton_raphson()
global COL ROW
global pw  sh  sw  sg  si cl  cm
global pcgw
global dnw  dng
global krw krg
global tsg tpcgw nrpermg
global eps
global T
global p_salinity T_salinity salinity
global INDC1 INDC2 INDC3
global  kx ky
global cnstx  dx dy dz
global k_temp  sgr swr wn gn  Pd0
global  lambda lambda_s lambda_w lambda_i lambda_h lambda_g phi
global pw_ini T_ini pw_scale T_scale cl_scale cm_scale


B=zeros(COL*ROW*4,1);
jac=spalloc(COL*ROW*4,COL*ROW*4,COL*ROW*4*20-16*8);

eps_pw=eps;
eps_cl=eps;
eps_T=eps;
lam_on=0;

pw_dimensionless=(pw-pw_ini)./pw_scale;
T_dimensionless=(T-T_ini)./T_scale;
cl_dimensionless=cl./cl_scale;
cm_dimensionless=cm./cm_scale;

for i=2:COL
    for j=1:ROW
        if INDC2(i,1)>1
            if INDC2(i,1)>3 && INDC3(i,1)==1
                cl(i,1)=interpolation2(p_salinity,T_salinity,salinity,pw(i,1)/1e6,T(i,1));
            elseif INDC3(i,1)>2
                cl(i,1)=((49.462^2-4*164.49*T(i,1))^0.5-49.462)/(2*164.49);
            elseif  INDC3(i,1)==2
                cl(i,1)=0;
            end
            if INDC3(i,1)==1 || (INDC3(i,1)>2 && INDC2(i,1)>3)
                if INDC1(i,1) == 0
                    cm(i,1)=methane_solubility(pw(i,1)/1e6,T(i,1),2*cl(i,1)/0.11099/(1-cl(i,1)))*0.016;
                else
                    cm(i,1)=hydrate_solubility(pw(i,1)/1e6,T(i,1),2*cl(i,1)/0.11099/(1-cl(i,1)))*0.016;
                end
            elseif INDC3(i,1)==2
                cm(i,1)=0;
            end

        end
        dnw(i,1)=brine_density(pw(i,1)/1e6,T(i,1),2*cl(i,1)/0.11099/(1-cl(i,1)),...
            cm(i,1)/0.016);
    end
end

for i=1:COL
    for j=1:ROW
        calc_trans(i,1);
    end
end




for i=1:COL

    % Define size of matrix.
    l=(i-1)*4+4;
    calc_trans(i,1);

    % Calculating residuals
    rsidw0=calc_Rw(i,1);
    rsidg0=calc_Rg(i,1);
    rsidc0=calc_Rc(i,1);
    rsidt0=calc_Rt(i,1);

    B(l-3)=-rsidw0;
    B(l-2)=-rsidg0;
    B(l-1)=-rsidc0;
    B(l)=-rsidt0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%          Code of primary variables switching                    %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%          1: pw, cm and cl                                       %%%%%
    %%%%%          2: pw, cm and cl -> pw, sw and cl                      %%%%%
    %%%%%          3: pw, sw and cl                                       %%%%%
    %%%%%          4: pw, sw and cl -> pw, sw and sh                      %%%%%
    %%%%%          5: pw, sw and sh                                       %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%                             i-1,1                              %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i > 1

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%          Change pw or Si         %%%%%%
        %%%%%%                                  %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if INDC2(i-1,1)<6

            pw(i-1,1)=(pw_dimensionless(i-1,1)+eps_pw)*pw_scale+pw_ini(i-1,1);
            temp1=dng(i-1,1);
            temp2=cl(i-1,1);
            temp3=cm(i-1,1);
            temp4=dnw(i-1,1);

            dng(i-1,1)=gas_density(pw(i-1,1)/1e6,T(i-1,1));
            if INDC2(i-1,1) > 3 && INDC3(i-1,1)==1
                cl(i-1,1)=interpolation2(p_salinity,T_salinity,salinity,pw(i-1,1)/1e6,T(i-1,1));
            end
            if (INDC2(i-1,1) > 1 && INDC3(i-1,1)==1)|| (INDC2(i-1,1)>3 && INDC3(i-1,1)>2)
                if INDC1(i-1,1) == 0
                    cm(i-1,1)=methane_solubility(pw(i-1,1)/1e6,T(i-1,1),2*cl(i-1,1)/0.11099/(1-cl(i-1,1)))*0.016;
                else
                    cm(i-1,1)=hydrate_solubility(pw(i-1,1)/1e6,T(i-1,1),2*cl(i-1,1)/0.11099/(1-cl(i-1,1)))*0.016;
                end
            end
            dnw(i-1,1)=brine_density(pw(i-1,1)/1e6,T(i-1,1),2*cl(i-1,1)/0.11099/(1-cl(i-1,1)),...
                cm(i-1,1)/0.016);
        else
            %*****************************************************************
            %********change pw to si as primary viarable**********************
            %*****************************************************************
            si(i-1,1)=si(i-1,1)+eps;
            temp1=lambda(i-1,1);
            temp2=kx(i-1,1);
            temp3=ky(i-1,1);
            temp4=cnstx(i-1,1);
            if i>2
                temp5=cnstx(i-2,1);
            end
            temp6=krg(i-1,1);
            temp7=krw(i-1,1);
            temp8=pcgw(i-1,1);

            if lam_on==1
                lambda(i-1,1)=(1-phi(i-1,1))*lambda_s+phi(i-1,1)*(sw(i-1,1)*lambda_w+(1-sw(i-1,1)-sh(i-1,1)-si(i-1,1))*lambda_g+sh(i-1,1)*lambda_h+si(i-1,1)*lambda_i);
            end

            if sh(i-1,1)>0 || si(i-1,1)>0
                tmp=(1-si(i-1,1)-sh(i-1,1))^2;
            else
                tmp=1;
            end
            kx(i-1,1)=k_temp(i-1,1)*tmp;
            ky(i-1,1)=k_temp(i-1,1)*tmp;

            cnstx(i-1,1)=2.0 * dz(i,1)*dy(i,1)*dz(i-1,1)*dy(i-1,1)*kx(i,1) * kx(i-1,1)/ ...
                (dx(i-1,1)*dz(i,1)*dy(i,1)*kx(i,1) + dx(i,1)*dz(i-1,1)*dy(i-1,1)*kx(i-1,1));
            if i>2
                cnstx(i-2,1)=2.0 * dz(i-1,1)*dy(i-1,1)*dz(i-2,1)*dy(i-2,1)*kx(i-1,1) * kx(i-2,1)/ ...
                    (dx(i-2,1)*dz(i-1,1)*dy(i-1,1)*kx(i-1,1) + dx(i-1,1)*dz(i-2,1)*dy(i-2,1)*kx(i-2,1));
            end

            if (1-sw(i-1,1)-sh(i-1,1)-si(i-1,1))/(1-sh(i-1,1)-si(i-1,1)) >= sgr
                krg(i-1,1)=((1-sw(i-1,1)-sh(i-1,1)-si(i-1,1))/(1-sh(i-1,1)-si(i-1,1))-sgr)^gn;
            else
                krg(i-1,1)=0;
            end
            if sw(i-1,1)/(1-sh(i-1,1)-si(i-1,1))>swr
                krw(i-1,1)=(sw(i-1,1)/(1-sh(i-1,1)-si(i-1,1))-swr)^wn;
            else
                krw(i-1,1)=0;
            end

            if sh(i-1,1)>0 || si(i-1,1)>0
                tmp_ratio=sqrt(1/(1+si(i-1,1)+sh(i-1,1)+2*(1-si(i-1,1)-sh(i-1,1))/log(si(i-1,1)+sh(i-1,1))));
            else
                tmp_ratio=1;
            end
            pcgw(i-1,1)=Pd0(i-1,1)*tmp_ratio*...
                interp(tsg,tpcgw,nrpermg,(1-sw(i-1,1)-sh(i-1,1)-si(i-1,1))/(1-sh(i-1,1)-si(i-1,1)),1);
        end

        calc_trans(i,1);
        rsidw1= calc_Rw(i,1);
        rsidg1= calc_Rg(i,1);
        rsidc1=calc_Rc(i,1);
        rsidt1=calc_Rt(i,1);
        jac(l-3,l-4*ROW-3)=(rsidw1-rsidw0)/eps_pw;
        jac(l-2,l-4*ROW-3)=(rsidg1-rsidg0)/eps_pw;
        jac(l-1,l-4*ROW-3)=(rsidc1-rsidc0)/eps_pw;
        jac(l,l-4*ROW-3)=(rsidt1-rsidt0)/eps_pw;

        if INDC2(i-1,1)<6
            pw(i-1,1)=pw_dimensionless(i-1,1)*pw_scale+pw_ini(i-1,1);
            dng(i-1,1)=temp1;
            cl(i-1,1)=temp2;
            cm(i-1,1)=temp3;
            dnw(i-1,1)=temp4;
        else
            %*****************************************************************
            %********change pw to si as primary viarable**********************
            %*****************************************************************
            si(i-1,1)=si(i-1,1)-eps;
            lambda(i-1,1)=temp1;
            kx(i-1,1)=temp2;
            ky(i-1,1)=temp3;
            cnstx(i-1,1)=temp4;
            if i>2
                cnstx(i-2,1)=temp5;
            end
            krg(i-1,1)=temp6;
            krw(i-1,1)=temp7;
            pcgw(i-1,1)=temp8;
        end
        calc_trans(i,1);



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%          Change sw     (or cm or si KH)          %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (INDC3(i-1,1)==1 && INDC2(i-1,1) == 1)
            cm(i-1,1)=(cm_dimensionless(i-1,1)+eps)*cm_scale;
            temp1=dnw(i-1,1);
            dnw(i-1,1)=brine_density(pw(i-1,1)/1e6,T(i-1,1),2*cl(i-1,1)/0.11099/(1-cl(i-1,1)),...
                cm(i-1,1)/0.016);
        elseif INDC3(i-1,1)~=2

            %*****************************************************************
            %********change cm to sw as primary viarable**********************
            %*****************************************************************

            sw(i-1,1)=sw(i-1,1)+eps;
            temp1=lambda(i-1,1);
            temp2=krw(i-1,1);
            temp3=krg(i-1,1);
            temp4=pcgw(i-1,1);

            if lam_on==1
                if INDC1(i-1,1)==0
                    lambda(i-1,1)=(1-phi(i-1,1))*lambda_s+phi(i-1,1)*(sw(i-1,1)*lambda_w+(1-sw(i-1,1)-sh(i-1,1)-si(i-1,1))*lambda_g+sh(i-1,1)*lambda_h+si(i-1,1)*lambda_i);
                else
                    lambda(i-1,1)=(1-phi(i-1,1))*lambda_s+phi(i-1,1)*(sw(i-1,1)*lambda_w+sh(i-1,1)*lambda_h+(1-sw(i-1,1)-sh(i-1,1))*lambda_i);
                end
            end

            if sw(i-1,1)/(1-sh(i-1,1)-si(i-1,1))>swr
                krw(i-1,1)=(sw(i-1,1)/(1-sh(i-1,1)-si(i-1,1))-swr)^wn;
            else
                krw(i-1,1)=0;
            end

            if INDC1(i-1,1)==1
                krg(i-1,1)=0;
            else
                if (1-sw(i-1,1)-sh(i-1,1)-si(i-1,1))/(1-sh(i-1,1)-si(i-1,1)) >= sgr
                    krg(i-1,1)=((1-sw(i-1,1)-sh(i-1,1)-si(i-1,1))/(1-sh(i-1,1)-si(i-1,1))-sgr)^gn;
                else
                    krg(i-1,1)=0;
                end
            end

            if sh(i-1,1)>0 || si(i-1,1)>0
                tmp_ratio=sqrt(1/(1+si(i-1,1)+sh(i-1,1)+2*(1-si(i-1,1)-sh(i-1,1))/log(si(i-1,1)+sh(i-1,1))));
            else
                tmp_ratio=1;
            end
            pcgw(i-1,1)=Pd0(i-1,1)*tmp_ratio*...
                interp(tsg,tpcgw,nrpermg,(1-sw(i-1,1)-sh(i-1,1)-si(i-1,1))/(1-sh(i-1,1)-si(i-1,1)),1);

        elseif INDC3(i-1,1)==2

            %*****************************************************************
            %********change sw to si as primary viarable**********************
            %*****************************************************************

            si(i-1,1)=si(i-1,1)+eps;
            temp1=lambda(i-1,1);
            temp2=kx(i-1,1);
            temp3=ky(i-1,1);
            temp4=cnstx(i-1,1);
            if i>2
                temp5=cnstx(i-2,1);
            end
            temp6=krg(i-1,1);
            temp7=krw(i-1,1);
            temp8=pcgw(i-1,1);


            if lam_on==1
                if INDC1(i-1,1)==1
                    lambda(i-1,1)=(1-phi(i-1,1))*lambda_s+phi(i-1,1)*((1-si(i-1,1))*lambda_h+si(i-1,1)*lambda_i);
                else
                    lambda(i-1,1)=(1-phi(i-1,1))*lambda_s+phi(i-1,1)*((1-si(i-1,1)-sh(i-1,1))*lambda_g+sh(i-1,1)*lambda_h+si(i-1,1)*lambda_i);
                end
            end

            if sh(i-1,1)>0 || si(i-1,1)>0
                tmp=(1-si(i-1,1)-sh(i-1,1))^2;
            else
                tmp=1;
            end
            kx(i-1,1)=k_temp(i-1,1)*tmp;
            ky(i-1,1)=k_temp(i-1,1)*tmp;

            cnstx(i-1,1)=2.0 * dz(i,1)*dy(i,1)*dz(i-1,1)*dy(i-1,1)*kx(i,1) * kx(i-1,1)/ ...
                (dx(i-1,1)*dz(i,1)*dy(i,1)*kx(i,1) + dx(i,1)*dz(i-1,1)*dy(i-1,1)*kx(i-1,1));
            if i>2
                cnstx(i-2,1)=2.0 * dz(i-1,1)*dy(i-1,1)*dz(i-2,1)*dy(i-2,1)*kx(i-1,1) * kx(i-2,1)/ ...
                    (dx(i-2,1)*dz(i-1,1)*dy(i-1,1)*kx(i-1,1) + dx(i-1,1)*dz(i-2,1)*dy(i-2,1)*kx(i-2,1));
            end

            if INDC1(i-1,1)==1
                krg(i-1,1)=0;
            else
                if (1-sh(i-1,1)-si(i-1,1))/(1-sh(i-1,1)-si(i-1,1)) >= sgr
                    krg(i-1,1)=((1-sh(i-1,1)-si(i-1,1))/(1-sh(i-1,1)-si(i-1,1))-sgr)^gn;
                else
                    krg(i-1,1)=0;
                end
            end
            krw(i-1,1)=0;

            if sh(i-1,1)>0 || si(i-1,1)>0
                tmp_ratio=sqrt(1/(1+si(i-1,1)+sh(i-1,1)+2*(1-si(i-1,1)-sh(i-1,1))/log(si(i-1,1)+sh(i-1,1))));
            else
                tmp_ratio=1;
            end
            pcgw(i-1,1)=Pd0(i-1,1)*tmp_ratio*...
                interp(tsg,tpcgw,nrpermg,(1-sh(i-1,1)-si(i-1,1))/(1-sh(i-1,1)-si(i-1,1)),1);
        end


        calc_trans(i,1);
        rsidw1= calc_Rw(i,1);
        rsidg1= calc_Rg(i,1);
        rsidc1=calc_Rc(i,1);
        rsidt1=calc_Rt(i,1);
        jac(l-3,l-4*ROW-2)=(rsidw1-rsidw0)/eps;
        jac(l-2,l-4*ROW-2)=(rsidg1-rsidg0)/eps;
        jac(l-1,l-4*ROW-2)=(rsidc1-rsidc0)/eps;
        jac(l,l-4*ROW-2)=(rsidt1-rsidt0)/eps;

        if (INDC3(i-1,1)==1 && INDC2(i-1,1) == 1)
            cm(i-1,1)=cm_dimensionless(i-1,1)*cm_scale;
            dnw(i-1,1)=temp1;
        elseif INDC3(i-1,1)~=2

            %*****************************************************************
            %********change cm to sw as primary viarable**********************
            %*****************************************************************

            sw(i-1,1)=sw(i-1,1)-eps;
            lambda(i-1,1)=temp1;
            krw(i-1,1)=temp2;
            krg(i-1,1)=temp3;
            pcgw(i-1,1)=temp4;

        elseif INDC3(i-1,1)==2

            %*****************************************************************
            %********change sw to si as primary viarable**********************
            %*****************************************************************

            si(i-1,1)=si(i-1,1)-eps;
            lambda(i-1,1)=temp1;
            kx(i-1,1)=temp2;
            ky(i-1,1)=temp3;
            cnstx(i-1,1)=temp4;
            if i>2
                cnstx(i-2,1)=temp5;
            end
            krg(i-1,1)=temp6;
            krw(i-1,1)=temp7;
            pcgw(i-1,1)=temp8;
        end

        calc_trans(i,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%          Change cl or sh or cm or si         %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if INDC2(i-1,1) <= 3 && INDC3(i-1,1)<3
            cl(i-1,1)=(cl_dimensionless(i-1,1)+eps)*cl_scale;
            temp1=dnw(i-1,1);

            dnw(i-1,1)=brine_density(pw(i-1,1)/1e6,T(i-1,1),2*cl(i-1,1)/0.11099/(1-cl(i-1,1)),...
                cm(i-1,1)/0.016);
        elseif INDC3(i-1,1)>2 && INDC2(i-1,1)<4

            %*****************************************************************
            %********change cl to cm as primary viarable**********************
            %*****************************************************************

            cm(i-1,1)=(cm_dimensionless(i-1,1)+eps)*cm_scale;
            temp1=dnw(i-1,1);

            dnw(i-1,1)=brine_density(pw(i-1,1)/1e6,T(i-1,1),2*cl(i-1,1)/0.11099/(1-cl(i-1,1)),...
                cm(i-1,1)/0.016);

        elseif INDC3(i-1,1)>2 && INDC1(i-1,1)==0 && INDC2(i-1,1)<6

            %*****************************************************************
            %********change cm to si as primary viarable**********************
            %*****************************************************************

            si(i-1,1)=si(i-1,1)+eps;
            temp1=lambda(i-1,1);
            temp2=kx(i-1,1);
            temp3=ky(i-1,1);
            temp4=cnstx(i-1,1);
            if i>2
                temp5=cnstx(i-2,1);
            end
            temp6=krg(i-1,1);
            temp7=krw(i-1,1);
            temp8=pcgw(i-1,1);

            if lam_on==1
                lambda(i-1,1)=(1-phi(i-1,1))*lambda_s+phi(i-1,1)*(sw(i-1,1)*lambda_w+(1-sw(i-1,1)-si(i-1,1))*lambda_g+si(i-1,1)*lambda_i);
            end

            if sh(i-1,1)>0 || si(i-1,1)>0
                tmp=(1-si(i-1,1)-sh(i-1,1))^2;
            else
                tmp=1;
            end
            kx(i-1,1)=k_temp(i-1,1)*tmp;
            ky(i-1,1)=k_temp(i-1,1)*tmp;

            cnstx(i-1,1)=2.0 * dz(i,1)*dy(i,1)*dz(i-1,1)*dy(i-1,1)*kx(i,1) * kx(i-1,1)/ ...
                (dx(i-1,1)*dz(i,1)*dy(i,1)*kx(i,1) + dx(i,1)*dz(i-1,1)*dy(i-1,1)*kx(i-1,1));
            if i>2
                cnstx(i-2,1)=2.0 * dz(i-1,1)*dy(i-1,1)*dz(i-2,1)*dy(i-2,1)*kx(i-1,1) * kx(i-2,1)/ ...
                    (dx(i-2,1)*dz(i-1,1)*dy(i-1,1)*kx(i-1,1) + dx(i-1,1)*dz(i-2,1)*dy(i-2,1)*kx(i-2,1));
            end


            if (1-sw(i-1,1)-si(i-1,1))/(1-sh(i-1,1)-si(i-1,1)) >= sgr
                krg(i-1,1)=((1-sw(i-1,1)-si(i-1,1))/(1-sh(i-1,1)-si(i-1,1))-sgr)^gn;
            else
                krg(i-1,1)=0;
            end

            if sw(i-1,1)/(1-sh(i-1,1)-si(i-1,1))>swr
                krw(i-1,1)=(sw(i-1,1)/(1-sh(i-1,1)-si(i-1,1))-swr)^wn;
            else
                krw(i-1,1)=0;
            end

            if sh(i-1,1)>0 || si(i-1,1)>0
                tmp_ratio=sqrt(1/(1+si(i-1,1)+sh(i-1,1)+2*(1-si(i-1,1)-sh(i-1,1))/log(si(i-1,1)+sh(i-1,1))));
            else
                tmp_ratio=1;
            end
            pcgw(i-1,1)=Pd0(i-1,1)*tmp_ratio*...
                interp(tsg,tpcgw,nrpermg,(1-sw(i-1,1)-si(i-1,1))/(1-sh(i-1,1)-si(i-1,1)),1);


        else

            %*****************************************************************
            %********change si to sh as primary viarable**********************
            %*****************************************************************

            sh(i-1,1)=sh(i-1,1)+eps;
            temp1=lambda(i-1,1);
            temp2=kx(i-1,1);
            temp3=ky(i-1,1);
            temp4=cnstx(i-1,1);
            if i>2
                temp5=cnstx(i-2,1);
            end
            temp6=krg(i-1,1);
            temp7=krw(i-1,1);
            temp8=pcgw(i-1,1);


            if lam_on==1
                if INDC1(i-1,1)==1
                    lambda(i-1,1)=(1-phi(i-1,1))*lambda_s+phi(i-1,1)*(sw(i-1,1)*lambda_w+sh(i-1,1)*lambda_h+(1-sh(i-1,1)-sw(i-1,1))*lambda_i);
                else
                    lambda(i-1,1)=(1-phi(i-1,1))*lambda_s+phi(i-1,1)*(sw(i-1,1)*lambda_w+(1-sw(i-1,1)-sh(i-1,1)-si(i-1,1))*lambda_g+sh(i-1,1)*lambda_h+si(i-1,1)*lambda_i);
                end
            end
            if sh(i-1,1)>0 || si(i-1,1)>0
                tmp=(1-si(i-1,1)-sh(i-1,1))^2;
            else
                tmp=1;
            end
            kx(i-1,1)=k_temp(i-1,1)*tmp;
            ky(i-1,1)=k_temp(i-1,1)*tmp;

            cnstx(i-1,1)=2.0 * dz(i,1)*dy(i,1)*dz(i-1,1)*dy(i-1,1)*kx(i,1) * kx(i-1,1)/ ...
                (dx(i-1,1)*dz(i,1)*dy(i,1)*kx(i,1) + dx(i,1)*dz(i-1,1)*dy(i-1,1)*kx(i-1,1));
            if i>2
                cnstx(i-2,1)=2.0 * dz(i-1,1)*dy(i-1,1)*dz(i-2,1)*dy(i-2,1)*kx(i-1,1) * kx(i-2,1)/ ...
                    (dx(i-2,1)*dz(i-1,1)*dy(i-1,1)*kx(i-1,1) + dx(i-1,1)*dz(i-2,1)*dy(i-2,1)*kx(i-2,1));
            end


            if sw(i-1,1)/(1-sh(i-1,1)-si(i-1,1))>swr
                krw(i-1,1)=(sw(i-1,1)/(1-sh(i-1,1)-si(i-1,1))-swr)^wn;
            else
                krw(i-1,1)=0;
            end

            if INDC1(i-1,1)==1
                krg(i-1,1)=0;
            else
                if (1-sw(i-1,1)-sh(i-1,1)-si(i-1,1))/(1-sh(i-1,1)-si(i-1,1)) >= sgr
                    krg(i-1,1)=((1-sw(i-1,1)-sh(i-1,1)-si(i-1,1))/(1-sh(i-1,1)-si(i-1,1))-sgr)^gn;
                else
                    krg(i-1,1)=0;
                end
            end

            if sh(i-1,1)>0 || si(i-1,1)>0
                tmp_ratio=sqrt(1/(1+si(i-1,1)+sh(i-1,1)+2*(1-si(i-1,1)-sh(i-1,1))/log(si(i-1,1)+sh(i-1,1))));
            else
                tmp_ratio=1;
            end
            pcgw(i-1,1)=Pd0(i-1,1)*tmp_ratio*...
                interp(tsg,tpcgw,nrpermg,(1-sw(i-1,1)-sh(i-1,1)-si(i-1,1))/(1-sh(i-1,1)-si(i-1,1)),1);

        end

        calc_trans(i,1);
        rsidw1= calc_Rw(i,1);
        rsidg1= calc_Rg(i,1);
        rsidc1=calc_Rc(i,1);
        rsidt1=calc_Rt(i,1);
        jac(l-3,l-4*ROW-1)=(rsidw1-rsidw0)/eps;
        jac(l-2,l-4*ROW-1)=(rsidg1-rsidg0)/eps;
        jac(l-1,l-4*ROW-1)=(rsidc1-rsidc0)/eps;
        jac(l,l-4*ROW-1)=(rsidt1-rsidt0)/eps;

        if INDC2(i-1,1) <= 3 && INDC3(i-1,1)<3
            cl(i-1,1)=cl_dimensionless(i-1,1)*cl_scale;
            dnw(i-1,1)=temp1;
        elseif INDC3(i-1,1)>2 && INDC2(i-1,1)<4

            %*****************************************************************
            %********change cl to cm as primary viarable**********************
            %*****************************************************************

            cm(i-1,1)=cm_dimensionless(i-1,1)*cm_scale;
            dnw(i-1,1)=temp1;

        elseif INDC3(i-1,1)>2 && INDC1(i-1,1)==0 && INDC2(i-1,1)<6

            %*****************************************************************
            %********change cm to si as primary viarable**********************
            %*****************************************************************

            si(i-1,1)=si(i-1,1)-eps;
            lambda(i-1,1)=temp1;
            kx(i-1,1)=temp2;
            ky(i-1,1)=temp3;
            cnstx(i-1,1)=temp4;
            if i>2
                cnstx(i-2,1)=temp5;
            end
            krg(i-1,1)=temp6;
            krw(i-1,1)=temp7;
            pcgw(i-1,1)=temp8;

        else

            %*****************************************************************
            %********change si to sh as primary viarable**********************
            %*****************************************************************

            sh(i-1,1)=sh(i-1,1)-eps;
            lambda(i-1,1)=temp1;
            kx(i-1,1)=temp2;
            ky(i-1,1)=temp3;
            cnstx(i-1,1)=temp4;
            if i>2
                cnstx(i-2,1)=temp5;
            end
            krg(i-1,1)=temp6;
            krw(i-1,1)=temp7;
            pcgw(i-1,1)=temp8;
        end

        calc_trans(i,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%          Change T                %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T(i-1,1)=(T_dimensionless(i-1,1)+eps_T)*T_scale+T_ini(i-1,1);
        temp1=cl(i-1,1);
        temp2=cm(i-1,1);
        temp3=dnw(i-1,1);
        temp4=dng(i-1,1);

        if INDC2(i-1,1) > 3 && INDC3(i-1,1)==1
            cl(i-1,1)=interpolation2(p_salinity,T_salinity,salinity,pw(i-1,1)/1e6,T(i-1,1));
        end
        if INDC3(i-1,1)>2
            cl(i-1,1)=((49.462^2-4*164.49*T(i-1,1))^0.5-49.462)/(2*164.49);
        end

        if (INDC2(i-1,1) > 1 && INDC3(i-1,1)==1) || (INDC3(i-1,1)>2 && INDC2(i-1,1)>3)
            if INDC1(i-1,1) == 0
                cm(i-1,1)=methane_solubility(pw(i-1,1)/1e6,T(i-1,1),2*cl(i-1,1)/0.11099/(1-cl(i-1,1)))*0.016;
            else
                cm(i-1,1)=hydrate_solubility(pw(i-1,1)/1e6,T(i-1,1),2*cl(i-1,1)/0.11099/(1-cl(i-1,1)))*0.016;
            end
        end

        dnw(i-1,1)=brine_density(pw(i-1,1)/1e6,T(i-1,1),2*cl(i-1,1)/0.11099/(1-cl(i-1,1)),...
            cm(i-1,1)/0.016);
        dng(i-1,1)=gas_density(pw(i-1,1)/1e6,T(i-1,1));

        calc_trans(i,1);
        rsidw1= calc_Rw(i,1);
        rsidg1= calc_Rg(i,1);
        rsidc1=calc_Rc(i,1);
        rsidt1=calc_Rt(i,1);
        jac(l-3,l-4*ROW)=(rsidw1-rsidw0)/eps_T;
        jac(l-2,l-4*ROW)=(rsidg1-rsidg0)/eps_T;
        jac(l-1,l-4*ROW)=(rsidc1-rsidc0)/eps_T;
        jac(l,l-4*ROW)=(rsidt1-rsidt0)/eps_T;

        T(i-1,1)=T_dimensionless(i-1,1)*T_scale+T_ini(i-1,1);
        cl(i-1,1)=temp1;
        cm(i-1,1)=temp2;
        dnw(i-1,1)=temp3;
        dng(i-1,1)=temp4;

        calc_trans(i,1);
    end % end of if i > 1



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%                              i, j                             %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%          Change pw or Sg         %%%%%%
    %%%%%%                                  %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if INDC2(i,1)<6
        pw(i,1)=(pw_dimensionless(i,1)+eps_pw)*pw_scale+pw_ini(i,1);
        temp1=dng(i,1);
        temp2=cl(i,1);
        temp3=cm(i,1);
        temp4=dnw(i,1);

        dng(i,1)=gas_density(pw(i,1)/1e6,T(i,1));
        if INDC2(i,1) > 3 && INDC3(i,1)==1
            cl(i,1)=interpolation2(p_salinity,T_salinity,salinity,pw(i,1)/1e6,T(i,1));
        end
        if (INDC2(i,1) > 1 && INDC3(i,1)==1)|| (INDC2(i,1)>3 && INDC3(i,1)>2)
            if INDC1(i,1) == 0
                cm(i,1)=methane_solubility(pw(i,1)/1e6,T(i,1),2*cl(i,1)/0.11099/(1-cl(i,1)))*0.016;
            else
                cm(i,1)=hydrate_solubility(pw(i,1)/1e6,T(i,1),2*cl(i,1)/0.11099/(1-cl(i,1)))*0.016;
            end
        end
        dnw(i,1)=brine_density(pw(i,1)/1e6,T(i,1),2*cl(i,1)/0.11099/(1-cl(i,1)),...
            cm(i,1)/0.016);
    else
        %*****************************************************************
        %********change pw to si as primary viarable**********************
        %*****************************************************************
        si(i,1)=si(i,1)+eps;
        temp1=lambda(i,1);
        temp2=kx(i,1);
        temp3=ky(i,1);
        if i<COL
            temp4=cnstx(i,1);
        end
        if i>1
            temp5=cnstx(i-1,1);
        end
        temp6=krg(i,1);
        temp7=krw(i,1);
        temp8=pcgw(i,1);

        if lam_on==1
            lambda(i,1)=(1-phi(i,1))*lambda_s+phi(i,1)*(sw(i,1)*lambda_w+(1-sw(i,1)-sh(i,1)-si(i,1))*lambda_g+sh(i,1)*lambda_h+si(i,1)*lambda_i);
        end

        if sh(i,1)>0 || si(i,1)>0
            tmp=(1-si(i,1)-sh(i,1))^2;
        else
            tmp=1;
        end
        kx(i,1)=k_temp(i,1)*tmp;
        ky(i,1)=k_temp(i,1)*tmp;

        if i < COL
            cnstx(i,1)=2.0 * dz(i+1,1)*dy(i+1,1)*dz(i,1)*dy(i,1)*kx(i+1,1) * kx(i,1)/ ...
                (dx(i,1)*dz(i+1,1)*dy(i+1,1)*kx(i+1,1) + dx(i+1,1)*dz(i,1)*dy(i,1)*kx(i,1));
        end
        if i>1
            cnstx(i-1,1)=2.0 * dz(i,1)*dy(i,1)*dz(i-1,1)*dy(i-1,1)*kx(i,1) * kx(i-1,1)/ ...
                (dx(i-1,1)*dz(i,1)*dy(i,1)*kx(i,1) + dx(i,1)*dz(i-1,1)*dy(i-1,1)*kx(i-1,1));
        end

        if (1-sw(i,1)-sh(i,1)-si(i,1))/(1-sh(i,1)-si(i,1)) >= sgr
            krg(i,1)=((1-sw(i,1)-sh(i,1)-si(i,1))/(1-sh(i,1)-si(i,1))-sgr)^gn;
        else
            krg(i,1)=0;
        end
        if sw(i,1)/(1-sh(i,1)-si(i,1))>swr
            krw(i,1)=(sw(i,1)/(1-sh(i,1)-si(i,1))-swr)^wn;
        else
            krw(i,1)=0;
        end

        if sh(i,1)>0 || si(i,1)>0
            tmp_ratio=sqrt(1/(1+si(i,1)+sh(i,1)+2*(1-si(i,1)-sh(i,1))/log(si(i,1)+sh(i,1))));
        else
            tmp_ratio=1;
        end
        pcgw(i,1)=Pd0(i,1)*tmp_ratio*...
            interp(tsg,tpcgw,nrpermg,(1-sw(i,1)-sh(i,1)-si(i,1))/(1-sh(i,1)-si(i,1)),1);
    end

    calc_trans(i,1);
    rsidw1= calc_Rw(i,1);
    rsidg1= calc_Rg(i,1);
    rsidc1=calc_Rc(i,1);
    rsidt1=calc_Rt(i,1);
    jac(l-3,l-3)=(rsidw1-rsidw0)/eps_pw;
    jac(l-2,l-3)=(rsidg1-rsidg0)/eps_pw;
    jac(l-1,l-3)=(rsidc1-rsidc0)/eps_pw;
    jac(l,l-3)=(rsidt1-rsidt0)/eps_pw;

    if INDC2(i,1)<6
        pw(i,1)=pw_dimensionless(i,1)*pw_scale+pw_ini(i,1);
        dng(i,1)=temp1;
        cl(i,1)=temp2;
        cm(i,1)=temp3;
        dnw(i,1)=temp4;
    else
        %*****************************************************************
        %********change pw to si as primary viarable**********************
        %*****************************************************************
        si(i,1)=si(i,1)-eps;
        lambda(i,1)=temp1;
        kx(i,1)=temp2;
        ky(i,1)=temp3;
        if i<COL
            cnstx(i,1)=temp4;
        end
        if i>1
            cnstx(i-1,1)=temp5;
        end
        krg(i,1)=temp6;
        krw(i,1)=temp7;
        pcgw(i,1)=temp8;
    end
    calc_trans(i,1);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%          Change sw     (or cm or si KH)          %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (INDC3(i,1)==1 && INDC2(i,1) == 1)
        cm(i,1)=(cm_dimensionless(i,1)+eps)*cm_scale;
        temp1=dnw(i,1);
        dnw(i,1)=brine_density(pw(i,1)/1e6,T(i,1),2*cl(i,1)/0.11099/(1-cl(i,1)),...
            cm(i,1)/0.016);
    elseif INDC3(i,1)~=2

        %*****************************************************************
        %********change cm to sw as primary viarable**********************
        %*****************************************************************

        sw(i,1)=sw(i,1)+eps;
        temp1=lambda(i,1);
        temp2=krw(i,1);
        temp3=krg(i,1);
        temp4=pcgw(i,1);

        if lam_on==1
            if INDC1(i,1)==0
                lambda(i,1)=(1-phi(i,1))*lambda_s+phi(i,1)*(sw(i,1)*lambda_w(i,1)+(1-sw(i,1)-sh(i,1)-si(i,1))*lambda_g+sh(i,1)*lambda_h+si(i,1)*lambda_i);
            else
                lambda(i,1)=(1-phi(i,1))*lambda_s+phi(i,1)*(sw(i,1)*lambda_w(i,1)+sh(i,1)*lambda_h+(1-sw(i,1)-sh(i,1))*lambda_i);
            end
        end

        if sw(i,1)/(1-sh(i,1)-si(i,1))>swr
            krw(i,1)=(sw(i,1)/(1-sh(i,1)-si(i,1))-swr)^wn;
        else
            krw(i,1)=0;
        end

        if INDC1(i,1)==1
            krg(i,1)=0;
        else
            if (1-sw(i,1)-sh(i,1)-si(i,1))/(1-sh(i,1)-si(i,1)) >= sgr
                krg(i,1)=((1-sw(i,1)-sh(i,1)-si(i,1))/(1-sh(i,1)-si(i,1))-sgr)^gn;
            else
                krg(i,1)=0;
            end
        end

        if sh(i,1)>0 || si(i,1)>0
            tmp_ratio=sqrt(1/(1+si(i,1)+sh(i,1)+2*(1-si(i,1)-sh(i,1))/log(si(i,1)+sh(i,1))));
        else
            tmp_ratio=1;
        end
        pcgw(i,1)=Pd0(i,1)*tmp_ratio*...
            interp(tsg,tpcgw,nrpermg,(1-sw(i,1)-sh(i,1)-si(i,1))/(1-sh(i,1)-si(i,1)),1);

    elseif INDC3(i,1)==2

        %*****************************************************************
        %********change sw to si as primary viarable**********************
        %*****************************************************************

        si(i,1)=si(i,1)+eps;
        temp1=lambda(i,1);
        temp2=kx(i,1);
        temp3=ky(i,1);
        if i<COL
            temp4=cnstx(i,1);
        end
        if i>1
            temp5=cnstx(i-1,1);
        end
        temp6=krg(i,1);
        temp7=krw(i,1);
        temp8=pcgw(i,1);


        if lam_on==1
            if INDC1(i,1)==1
                lambda(i,1)=(1-phi(i,1))*lambda_s+phi(i,1)*((1-si(i,1))*lambda_h+si(i,1)*lambda_i);
            else
                lambda(i,1)=(1-phi(i,1))*lambda_s+phi(i,1)*((1-si(i,1)-sh(i,1))*lambda_g+sh(i,1)*lambda_h+si(i,1)*lambda_i);
            end
        end

        if sh(i,1)>0 || si(i,1)>0
            tmp=(1-si(i,1)-sh(i,1))^2;
        else
            tmp=1;
        end
        kx(i,1)=k_temp(i,1)*tmp;
        ky(i,1)=k_temp(i,1)*tmp;

        if i < COL
            cnstx(i,1)=2.0 * dz(i+1,1)*dy(i+1,1)*dz(i,1)*dy(i,1)*kx(i+1,1) * kx(i,1)/ ...
                (dx(i,1)*dz(i+1,1)*dy(i+1,1)*kx(i+1,1) + dx(i+1,1)*dz(i,1)*dy(i,1)*kx(i,1));
        end
        if i>1
            cnstx(i-1,1)=2.0 * dz(i,1)*dy(i,1)*dz(i-1,1)*dy(i-1,1)*kx(i,1) * kx(i-1,1)/ ...
                (dx(i-1,1)*dz(i,1)*dy(i,1)*kx(i,1) + dx(i,1)*dz(i-1,1)*dy(i-1,1)*kx(i-1,1));
        end

        if INDC1(i,1)==1
            krg(i,1)=0;
        else
            if (1-sh(i,1)-si(i,1))/(1-sh(i,1)-si(i,1)) >= sgr
                krg(i,1)=((1-sh(i,1)-si(i,1))/(1-sh(i,1)-si(i,1))-sgr)^gn;
            else
                krg(i,1)=0;
            end
        end
        krw(i,1)=0;

        if sh(i,1)>0 || si(i,1)>0
            tmp_ratio=sqrt(1/(1+si(i,1)+sh(i,1)+2*(1-si(i,1)-sh(i,1))/log(si(i,1)+sh(i,1))));
        else
            tmp_ratio=1;
        end
        pcgw(i,1)=Pd0(i,1)*tmp_ratio*...
            interp(tsg,tpcgw,nrpermg,(1-sh(i,1)-si(i,1))/(1-sh(i,1)-si(i,1)),1);
    end


    calc_trans(i,1);
    rsidw1= calc_Rw(i,1);
    rsidg1= calc_Rg(i,1);
    rsidc1=calc_Rc(i,1);
    rsidt1=calc_Rt(i,1);
    jac(l-3,l-2)=(rsidw1-rsidw0)/eps;
    jac(l-2,l-2)=(rsidg1-rsidg0)/eps;
    jac(l-1,l-2)=(rsidc1-rsidc0)/eps;
    jac(l,l-2)=(rsidt1-rsidt0)/eps;

    if (INDC3(i,1)==1 && INDC2(i,1) == 1)
        cm(i,1)=cm_dimensionless(i,1)*cm_scale;
        dnw(i,1)=temp1;
    elseif INDC3(i,1)~=2

        %*****************************************************************
        %********change cm to sw as primary viarable**********************
        %*****************************************************************

        sw(i,1)=sw(i,1)-eps;
        lambda(i,1)=temp1;
        krw(i,1)=temp2;
        krg(i,1)=temp3;
        pcgw(i,1)=temp4;

    elseif INDC3(i,1)==2

        %*****************************************************************
        %********change sw to si as primary viarable**********************
        %*****************************************************************

        si(i,1)=si(i,1)-eps;
        lambda(i,1)=temp1;
        kx(i,1)=temp2;
        ky(i,1)=temp3;
        if i<COL
            cnstx(i,1)=temp4;
        end
        if i>1
            cnstx(i-1,1)=temp5;
        end
        krg(i,1)=temp6;
        krw(i,1)=temp7;
        pcgw(i,1)=temp8;
    end

    calc_trans(i,1);




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%          Change cl or sh or cm or si         %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if INDC2(i,1) <= 3 && INDC3(i,1)<3
        cl(i,1)=(cl_dimensionless(i,1)+eps)*cl_scale;
        temp1=dnw(i,1);
        dnw(i,1)=brine_density(pw(i,1)/1e6,T(i,1),2*cl(i,1)/0.11099/(1-cl(i,1)),...
            cm(i,1)/0.016);
    elseif INDC3(i,1)>2 && INDC2(i,1)<4

        %*****************************************************************
        %********change cl to cm as primary viarable**********************
        %*****************************************************************

        cm(i,1)=(cm_dimensionless(i,1)+eps)*cm_scale;
        temp1=dnw(i,1);
        dnw(i,1)=brine_density(pw(i,1)/1e6,T(i,1),2*cl(i,1)/0.11099/(1-cl(i,1)),...
            cm(i,1)/0.016);

    elseif INDC3(i,1)>2 && INDC1(i,1)==0 && INDC2(i,1)<6

        %*****************************************************************
        %********change cm to si as primary viarable**********************
        %*****************************************************************

        si(i,1)=si(i,1)+eps;
        temp1=lambda(i,1);
        temp2=kx(i,1);
        temp3=ky(i,1);
        if i<COL
            temp4=cnstx(i,1);
        end
        if i>1
            temp5=cnstx(i-1,1);
        end
        temp6=krg(i,1);
        temp7=krw(i,1);
        temp8=pcgw(i,1);


        if lam_on==1
            lambda(i,1)=(1-phi(i,1))*lambda_s+phi(i,1)*(sw(i,1)*lambda_w(i,1)+(1-sw(i,1)-si(i,1))*lambda_g+sh(i,1)*lambda_h+si(i,1)*lambda_i);
        end

        if sh(i,1)>0 || si(i,1)>0
            tmp=(1-si(i,1)-sh(i,1))^2;
        else
            tmp=1;
        end
        kx(i,1)=k_temp(i,1)*tmp;
        ky(i,1)=k_temp(i,1)*tmp;

        if i < COL
            cnstx(i,1)=2.0 * dz(i+1,1)*dy(i+1,1)*dz(i,1)*dy(i,1)*kx(i+1,1) * kx(i,1)/ ...
                (dx(i,1)*dz(i+1,1)*dy(i+1,1)*kx(i+1,1) + dx(i+1,1)*dz(i,1)*dy(i,1)*kx(i,1));
        end
        if i>1
            cnstx(i-1,1)=2.0 * dz(i,1)*dy(i,1)*dz(i-1,1)*dy(i-1,1)*kx(i,1) * kx(i-1,1)/ ...
                (dx(i-1,1)*dz(i,1)*dy(i,1)*kx(i,1) + dx(i,1)*dz(i-1,1)*dy(i-1,1)*kx(i-1,1));
        end


        if (1-sw(i,1)-si(i,1))/(1-sh(i,1)-si(i,1)) >= sgr
            krg(i,1)=((1-sw(i,1)-si(i,1))/(1-sh(i,1)-si(i,1))-sgr)^gn;
        else
            krg(i,1)=0;
        end

        if sw(i,1)/(1-sh(i,1)-si(i,1))>swr
            krw(i,1)=(sw(i,1)/(1-sh(i,1)-si(i,1))-swr)^wn;
        else
            krw(i,1)=0;
        end

        if sh(i,1)>0 || si(i,1)>0
            tmp_ratio=sqrt(1/(1+si(i,1)+sh(i,1)+2*(1-si(i,1)-sh(i,1))/log(si(i,1)+sh(i,1))));
        else
            tmp_ratio=1;
        end
        pcgw(i,1)=Pd0(i,1)*tmp_ratio*...
            interp(tsg,tpcgw,nrpermg,(1-sw(i,1)-si(i,1))/(1-sh(i,1)-si(i,1)),1);


    else

        %*****************************************************************
        %********change si to sh as primary viarable**********************
        %*****************************************************************

        sh(i,1)=sh(i,1)+eps;
        temp1=lambda(i,1);
        temp2=kx(i,1);
        temp3=ky(i,1);
        if i<COL
            temp4=cnstx(i,1);
        end
        if i>1
            temp5=cnstx(i-1,1);
        end
        temp6=krg(i,1);
        temp7=krw(i,1);
        temp8=pcgw(i,1);

        if lam_on==1
            if INDC1(i,1)==1
                lambda(i,1)=(1-phi(i,1))*lambda_s+phi(i,1)*(sw(i,1)*lambda_w+sh(i,1)*lambda_h+(1-sh(i,1)-sw(i,1))*lambda_i);
            else
                lambda(i,1)=(1-phi(i,1))*lambda_s+phi(i,1)*(sw(i,1)*lambda_w+(1-sw(i,1)-sh(i,1)-si(i,1))*lambda_g+sh(i,1)*lambda_h+si(i,1)*lambda_i);
            end
        end
        if sh(i,1)>0 || si(i,1)>0
            tmp=(1-si(i,1)-sh(i,1))^2;
        else
            tmp=1;
        end
        kx(i,1)=k_temp(i,1)*tmp;
        ky(i,1)=k_temp(i,1)*tmp;

        if i < COL
            cnstx(i,1)=2.0 * dz(i+1,1)*dy(i+1,1)*dz(i,1)*dy(i,1)*kx(i+1,1) * kx(i,1)/ ...
                (dx(i,1)*dz(i+1,1)*dy(i+1,1)*kx(i+1,1) + dx(i+1,1)*dz(i,1)*dy(i,1)*kx(i,1));
        end
        if i>1
            cnstx(i-1,1)=2.0 * dz(i,1)*dy(i,1)*dz(i-1,1)*dy(i-1,1)*kx(i,1) * kx(i-1,1)/ ...
                (dx(i-1,1)*dz(i,1)*dy(i,1)*kx(i,1) + dx(i,1)*dz(i-1,1)*dy(i-1,1)*kx(i-1,1));
        end


        if sw(i,1)/(1-sh(i,1)-si(i,1))>swr
            krw(i,1)=(sw(i,1)/(1-sh(i,1)-si(i,1))-swr)^wn;
        else
            krw(i,1)=0;
        end

        if INDC1(i,1)==1
            krg(i,1)=0;
        else
            if (1-sw(i,1)-sh(i,1)-si(i,1))/(1-sh(i,1)-si(i,1)) >= sgr
                krg(i,1)=((1-sw(i,1)-sh(i,1)-si(i,1))/(1-sh(i,1)-si(i,1))-sgr)^gn;
            else
                krg(i,1)=0;
            end
        end

        if sh(i,1)>0 || si(i,1)>0
            tmp_ratio=sqrt(1/(1+si(i,1)+sh(i,1)+2*(1-si(i,1)-sh(i,1))/log(si(i,1)+sh(i,1))));
        else
            tmp_ratio=1;
        end
        pcgw(i,1)=Pd0(i,1)*tmp_ratio*...
            interp(tsg,tpcgw,nrpermg,(1-sw(i,1)-sh(i,1)-si(i,1))/(1-sh(i,1)-si(i,1)),1);
    end

    calc_trans(i,1);
    rsidw1= calc_Rw(i,1);
    rsidg1= calc_Rg(i,1);
    rsidc1=calc_Rc(i,1);
    rsidt1=calc_Rt(i,1);
    jac(l-3,l-1)=(rsidw1-rsidw0)/eps;
    jac(l-2,l-1)=(rsidg1-rsidg0)/eps;
    jac(l-1,l-1)=(rsidc1-rsidc0)/eps;
    jac(l,l-1)=(rsidt1-rsidt0)/eps;

    if INDC2(i,1) <= 3 && INDC3(i,1)<3
        cl(i,1)=cl_dimensionless(i,1)*cl_scale;
        dnw(i,1)=temp1;
    elseif INDC3(i,1)>2 && INDC2(i,1)<4

        %*****************************************************************
        %********change cl to cm as primary viarable**********************
        %*****************************************************************

        cm(i,1)=cm_dimensionless(i,1)*cm_scale;
        dnw(i,1)=temp1;

    elseif INDC3(i,1)>2 && INDC1(i,1)==0 && INDC2(i,1)<6

        %*****************************************************************
        %********change cm to si as primary viarable**********************
        %*****************************************************************

        si(i,1)=si(i,1)-eps;
        lambda(i,1)=temp1;
        kx(i,1)=temp2;
        ky(i,1)=temp3;
        if i<COL
            cnstx(i,1)=temp4;
        end
        if i>1
            cnstx(i-1,1)=temp5;
        end
        krg(i,1)=temp6;
        krw(i,1)=temp7;
        pcgw(i,1)=temp8;
    else

        %*****************************************************************
        %********change si to sh as primary viarable**********************
        %*****************************************************************

        sh(i,1)=sh(i,1)-eps;
        lambda(i,1)=temp1;
        kx(i,1)=temp2;
        ky(i,1)=temp3;
        if i<COL
            cnstx(i,1)=temp4;
        end
        if i>1
            cnstx(i-1,1)=temp5;
        end
        krg(i,1)=temp6;
        krw(i,1)=temp7;
        pcgw(i,1)=temp8;
    end
    calc_trans(i,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%          Change T                %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T(i,1)=(T_dimensionless(i,1)+eps_T)*T_scale+T_ini(i,1);
    temp1=cl(i,1);
    temp2=cm(i,1);
    temp3=dnw(i,1);
    temp4=dng(i,1);

    if INDC2(i,1) > 3 && INDC3(i,1)==1
        cl(i,1)=interpolation2(p_salinity,T_salinity,salinity,pw(i,1)/1e6,T(i,1));
    end
    if INDC3(i,1)>2
        cl(i,1)=((49.462^2-4*164.49*T(i,1))^0.5-49.462)/(2*164.49);
    end

    if (INDC2(i,1) > 1 && INDC3(i,1)==1) || (INDC3(i,1)>2 && INDC2(i,1)>3)
        if INDC1(i,1) == 0
            cm(i,1)=methane_solubility(pw(i,1)/1e6,T(i,1),2*cl(i,1)/0.11099/(1-cl(i,1)))*0.016;
        else
            cm(i,1)=hydrate_solubility(pw(i,1)/1e6,T(i,1),2*cl(i,1)/0.11099/(1-cl(i,1)))*0.016;
        end
    end

    dnw(i,1)=brine_density(pw(i,1)/1e6,T(i,1),2*cl(i,1)/0.11099/(1-cl(i,1)),...
        cm(i,1)/0.016);
    dng(i,1)=gas_density(pw(i,1)/1e6,T(i,1));

    calc_trans(i,1);
    rsidw1= calc_Rw(i,1);
    rsidg1= calc_Rg(i,1);
    rsidc1=calc_Rc(i,1);
    rsidt1=calc_Rt(i,1);
    jac(l-3,l)=(rsidw1-rsidw0)/eps_T;
    jac(l-2,l)=(rsidg1-rsidg0)/eps_T;
    jac(l-1,l)=(rsidc1-rsidc0)/eps_T;
    jac(l,l)=(rsidt1-rsidt0)/eps_T;

    T(i,1)=T_dimensionless(i,1)*T_scale+T_ini(i,1);
    cl(i,1)=temp1;
    cm(i,1)=temp2;
    dnw(i,1)=temp3;
    dng(i,1)=temp4;
    calc_trans(i,1);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%           i+1, j                %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i < COL
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%          Change pw or Si         %%%%%%
        %%%%%%                                  %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if INDC2(i+1,1)<6
            pw(i+1,1)=(pw_dimensionless(i+1,1)+eps_pw)*pw_scale+pw_ini(i+1,1);
            temp1=dng(i+1,1);
            temp2=cl(i+1,1);
            temp3=cm(i+1,1);
            temp4=dnw(i+1,1);

            dng(i+1,1)=gas_density(pw(i+1,1)/1e6,T(i+1,1));
            if INDC2(i+1,1) > 3 && INDC3(i+1,1)==1
                cl(i+1,1)=interpolation2(p_salinity,T_salinity,salinity,pw(i+1,1)/1e6,T(i+1,1));
            end
            if (INDC2(i+1,1) > 1 && INDC3(i+1,1)==1)|| (INDC2(i+1,1)>3 && INDC3(i+1,1)>2)
                if INDC1(i+1,1) == 0
                    cm(i+1,1)=methane_solubility(pw(i+1,1)/1e6,T(i+1,1),2*cl(i+1,1)/0.11099/(1-cl(i+1,1)))*0.016;
                else
                    cm(i+1,1)=hydrate_solubility(pw(i+1,1)/1e6,T(i+1,1),2*cl(i+1,1)/0.11099/(1-cl(i+1,1)))*0.016;
                end
            end
            dnw(i+1,1)=brine_density(pw(i+1,1)/1e6,T(i+1,1),2*cl(i+1,1)/0.11099/(1-cl(i+1,1)),...
                cm(i+1,1)/0.016);
        else
            %*****************************************************************
            %********change pw to si as primary viarable**********************
            %*****************************************************************
            si(i+1,1)=si(i+1,1)+eps;
            temp1=lambda(i+1,1);
            temp2=kx(i+1,1);
            temp3=ky(i+1,1);
            if i+1<COL
                temp4=cnstx(i+1,1);
            end

            if i<COL
                temp5=cnstx(i,1);
            end
            temp6=krg(i+1,1);
            temp7=krw(i+1,1);
            temp8=pcgw(i+1,1);

            if lam_on==1
                lambda(i+1,1)=(1-phi(i+1,1))*lambda_s+phi(i+1,1)*(sw(i+1,1)*lambda_w+(1-sw(i+1,1)-sh(i+1,1)-si(i+1,1))*lambda_g+sh(i+1,1)*lambda_h+si(i+1,1)*lambda_i);
            end

            if sh(i+1,1)>0 || si(i+1,1)>0
                tmp=(1-si(i+1,1)-sh(i+1,1))^2;
            else
                tmp=1;
            end
            kx(i+1,1)=k_temp(i+1,1)*tmp;
            ky(i+1,1)=k_temp(i+1,1)*tmp;

            if (i+1) < COL
                cnstx(i+1,1)=2.0 * dz(i+2,1)*dy(i+2,1)*dz(i+1,1)*dy(i+1,1)*kx(i+2,1) * kx(i+1,1)/ ...
                    (dx(i+1,1)*dz(i+2,1)*dy(i+2,1)*kx(i+2,1) + dx(i+2,1)*dz(i+1,1)*dy(i+1,1)*kx(i+1,1));
            end
            if i < COL
                cnstx(i,1)=2.0 * dz(i+1,1)*dy(i+1,1)*dz(i,1)*dy(i,1)*kx(i+1,1) * kx(i,1)/ ...
                    (dx(i,1)*dz(i+1,1)*dy(i+1,1)*kx(i+1,1) + dx(i+1,1)*dz(i,1)*dy(i,1)*kx(i,1));
            end

            if (1-sw(i+1,1)-sh(i+1,1)-si(i+1,1))/(1-sh(i+1,1)-si(i+1,1)) >= sgr
                krg(i+1,1)=((1-sw(i+1,1)-sh(i+1,1)-si(i+1,1))/(1-sh(i+1,1)-si(i+1,1))-sgr)^gn;
            else
                krg(i+1,1)=0;
            end
            if sw(i+1,1)/(1-sh(i+1,1)-si(i+1,1))>swr
                krw(i+1,1)=(sw(i+1,1)/(1-sh(i+1,1)-si(i+1,1))-swr)^wn;
            else
                krw(i+1,1)=0;
            end

            if sh(i+1,1)>0 || si(i+1,1)>0
                tmp_ratio=sqrt(1/(1+si(i+1,1)+sh(i+1,1)+2*(1-si(i+1,1)-sh(i+1,1))/log(si(i+1,1)+sh(i+1,1))));
            else
                tmp_ratio=1;
            end
            pcgw(i+1,1)=Pd0(i+1,1)*tmp_ratio*...
                interp(tsg,tpcgw,nrpermg,(1-sw(i+1,1)-sh(i+1,1)-si(i+1,1))/(1-sh(i+1,1)-si(i+1,1)),1);
        end

        calc_trans(i,1);
        rsidw1= calc_Rw(i,1);
        rsidg1= calc_Rg(i,1);
        rsidc1=calc_Rc(i,1);
        rsidt1=calc_Rt(i,1);
        jac(l-3,l+4*ROW-3)=(rsidw1-rsidw0)/eps_pw;
        jac(l-2,l+4*ROW-3)=(rsidg1-rsidg0)/eps_pw;
        jac(l-1,l+4*ROW-3)=(rsidc1-rsidc0)/eps_pw;
        jac(l,l+4*ROW-3)=(rsidt1-rsidt0)/eps_pw;

        if INDC2(i+1,1)<6
            pw(i+1,1)=pw_dimensionless(i+1,1)*pw_scale+pw_ini(i+1,1);
            dng(i+1,1)=temp1;
            cl(i+1,1)=temp2;
            cm(i+1,1)=temp3;
            dnw(i+1,1)=temp4;
        else
            %*****************************************************************
            %********change pw to si as primary viarable**********************
            %*****************************************************************
            si(i+1,1)=si(i+1,1)-eps;
            lambda(i+1,1)=temp1;
            kx(i+1,1)=temp2;
            ky(i+1,1)=temp3;
            if i+1<COL
                cnstx(i+1,1)=temp4;
            end
            if i<COL
                cnstx(i,1)=temp5;
            end
            krg(i+1,1)=temp6;
            krw(i+1,1)=temp7;
            pcgw(i+1,1)=temp8;
        end
        calc_trans(i,1);




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%          Change sw     (or cm or si KH)          %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (INDC3(i+1,1)==1 && INDC2(i+1,1) == 1)
            cm(i+1,1)=(cm_dimensionless(i+1,1)+eps)*cm_scale;
            temp1=dnw(i+1,1);
            dnw(i+1,1)=brine_density(pw(i+1,1)/1e6,T(i+1,1),2*cl(i+1,1)/0.11099/(1-cl(i+1,1)),...
                cm(i+1,1)/0.016);
        elseif INDC3(i+1,1)~=2

            %*****************************************************************
            %********change cm to sw as primary viarable**********************
            %*****************************************************************

            sw(i+1,1)=sw(i+1,1)+eps;
            temp1=lambda(i+1,1);
            temp2=krw(i+1,1);
            temp3=krg(i+1,1);
            temp4=pcgw(i+1,1);

            if lam_on==1
                if INDC1(i+1,1)==0
                    lambda(i+1,1)=(1-phi(i+1,1))*lambda_s+phi(i+1,1)*(sw(i+1,1)*lambda_w+(1-sw(i+1,1)-sh(i+1,1)-si(i+1,1))*lambda_g+sh(i+1,1)*lambda_h+si(i+1,1)*lambda_i);
                else
                    lambda(i+1,1)=(1-phi(i+1,1))*lambda_s+phi(i+1,1)*(sw(i+1,1)*lambda_w+(1-sw(i+1,1)-sh(i+1,1))*lambda_i+sh(i+1,1)*lambda_h);
                end
            end

            if sw(i+1,1)/(1-sh(i+1,1)-si(i+1,1))>swr
                krw(i+1,1)=(sw(i+1,1)/(1-sh(i+1,1)-si(i+1,1))-swr)^wn;
            else
                krw(i+1,1)=0;
            end

            if INDC1(i+1,1)==1
                krg(i+1,1)=0;
            else
                if (1-sw(i+1,1)-sh(i+1,1)-si(i+1,1))/(1-sh(i+1,1)-si(i+1,1)) >= sgr
                    krg(i+1,1)=((1-sw(i+1,1)-sh(i+1,1)-si(i+1,1))/(1-sh(i+1,1)-si(i+1,1))-sgr)^gn;
                else
                    krg(i+1,1)=0;
                end
            end

            if sh(i+1,1)>0 || si(i+1,1)>0
                tmp_ratio=sqrt(1/(1+si(i+1,1)+sh(i+1,1)+2*(1-si(i+1,1)-sh(i+1,1))/log(si(i+1,1)+sh(i+1,1))));
            else
                tmp_ratio=1;
            end
            pcgw(i+1,1)=Pd0(i+1,1)*tmp_ratio*...
                interp(tsg,tpcgw,nrpermg,(1-sw(i+1,1)-sh(i+1,1)-si(i+1,1))/(1-sh(i+1,1)-si(i+1,1)),1);

        elseif INDC3(i+1,1)==2

            %*****************************************************************
            %********change sw to si as primary viarable**********************
            %*****************************************************************

            si(i+1,1)=si(i+1,1)+eps;
            temp1=lambda(i+1,1);
            temp2=kx(i+1,1);
            temp3=ky(i+1,1);
            if i+1<COL
                temp4=cnstx(i+1,1);
            end
            if i<COL
                temp5=cnstx(i,1);
            end
            temp6=krg(i+1,1);
            temp7=krw(i+1,1);
            temp8=pcgw(i+1,1);

            if lam_on==1
                if INDC1(i+1,1)==1
                    lambda(i+1,1)=(1-phi(i+1,1))*lambda_s+phi(i+1,1)*((1-si(i+1,1))*lambda_h+si(i+1,1)*lambda_i);
                else
                    lambda(i+1,1)=(1-phi(i+1,1))*lambda_s+phi(i+1,1)*((1-si(i+1,1)-sh(i+1,1))*lambda_g+sh(i+1,1)*lambda_h+si(i+1,1)*lambda_i);
                end
            end

            if sh(i+1,1)>0 || si(i+1,1)>0
                tmp=(1-si(i+1,1)-sh(i+1,1))^2;
            else
                tmp=1;
            end
            kx(i+1,1)=k_temp(i+1,1)*tmp;
            ky(i+1,1)=k_temp(i+1,1)*tmp;

            if (i+1) < COL
                cnstx(i+1,1)=2.0 * dz(i+2,1)*dy(i+2,1)*dz(i+1,1)*dy(i+1,1)*kx(i+2,1) * kx(i+1,1)/ ...
                    (dx(i+1,1)*dz(i+2,1)*dy(i+2,1)*kx(i+2,1) + dx(i+2,1)*dz(i+1,1)*dy(i+1,1)*kx(i+1,1));
            end
            if i < COL
                cnstx(i,1)=2.0 * dz(i+1,1)*dy(i+1,1)*dz(i,1)*dy(i,1)*kx(i+1,1) * kx(i,1)/ ...
                    (dx(i,1)*dz(i+1,1)*dy(i+1,1)*kx(i+1,1) + dx(i+1,1)*dz(i,1)*dy(i,1)*kx(i,1));
            end

            if INDC1(i+1,1)==1
                krg(i+1,1)=0;
            else
                if (1-sh(i+1,1)-si(i+1,1))/(1-sh(i+1,1)-si(i+1,1)) >= sgr
                    krg(i+1,1)=((1-sh(i+1,1)-si(i+1,1))/(1-sh(i+1,1)-si(i+1,1))-sgr)^gn;
                else
                    krg(i+1,1)=0;
                end
            end
            krw(i+1,1)=0;

            if sh(i+1,1)>0 || si(i+1,1)>0
                tmp_ratio=sqrt(1/(1+si(i+1,1)+sh(i+1,1)+2*(1-si(i+1,1)-sh(i+1,1))/log(si(i+1,1)+sh(i+1,1))));
            else
                tmp_ratio=1;
            end
            pcgw(i+1,1)=Pd0(i+1,1)*tmp_ratio*...
                interp(tsg,tpcgw,nrpermg,(1-sh(i+1,1)-si(i+1,1))/(1-sh(i+1,1)-si(i+1,1)),1);
        end


        calc_trans(i,1);
        rsidw1= calc_Rw(i,1);
        rsidg1= calc_Rg(i,1);
        rsidc1=calc_Rc(i,1);
        rsidt1=calc_Rt(i,1);
        jac(l-3,l+4*ROW-2)=(rsidw1-rsidw0)/eps;
        jac(l-2,l+4*ROW-2)=(rsidg1-rsidg0)/eps;
        jac(l-1,l+4*ROW-2)=(rsidc1-rsidc0)/eps;
        jac(l,l+4*ROW-2)=(rsidt1-rsidt0)/eps;

        if (INDC3(i+1,1)==1 && INDC2(i+1,1) == 1)
            cm(i+1,1)=cm_dimensionless(i+1,1)*cm_scale;
            dnw(i+1,1)=temp1;
        elseif INDC3(i+1,1)~=2

            %*****************************************************************
            %********change cm to sw as primary viarable**********************
            %*****************************************************************

            sw(i+1,1)=sw(i+1,1)-eps;
            lambda(i+1,1)=temp1;
            krw(i+1,1)=temp2;
            krg(i+1,1)=temp3;
            pcgw(i+1,1)=temp4;

        elseif INDC3(i+1,1)==2

            %*****************************************************************
            %********change sw to si as primary viarable**********************
            %*****************************************************************

            si(i+1,1)=si(i+1,1)-eps;
            lambda(i+1,1)=temp1;
            kx(i+1,1)=temp2;
            ky(i+1,1)=temp3;
            if i+1<COL
                cnstx(i+1,1)=temp4;
            end
            if i<COL
                cnstx(i,1)=temp5;
            end
            krg(i+1,1)=temp6;
            krw(i+1,1)=temp7;
            pcgw(i+1,1)=temp8;
        end

        calc_trans(i,1);




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%          Change cl or sh or cm or si         %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if INDC2(i+1,1) <= 3 && INDC3(i+1,1)<3
            cl(i+1,1)=(cl_dimensionless(i+1,1)+eps)*cl_scale;
            temp1=dnw(i+1,1);
            dnw(i+1,1)=brine_density(pw(i+1,1)/1e6,T(i+1,1),2*cl(i+1,1)/0.11099/(1-cl(i+1,1)),...
                cm(i+1,1)/0.016);
        elseif INDC3(i+1,1)>2 && INDC2(i+1,1)<4

            %*****************************************************************
            %********change cl to cm as primary viarable**********************
            %*****************************************************************

            cm(i+1,1)=(cm_dimensionless(i+1,1)+eps)*cm_scale;
            temp1=dnw(i+1,1);
            dnw(i+1,1)=brine_density(pw(i+1,1)/1e6,T(i+1,1),2*cl(i+1,1)/0.11099/(1-cl(i+1,1)),...
                cm(i+1,1)/0.016);

        elseif INDC3(i+1,1)>2 && INDC1(i+1,1)==0 && INDC2(i+1,1)<6

            %*****************************************************************
            %********change cm to si as primary viarable**********************
            %*****************************************************************

            si(i+1,1)=si(i+1,1)+eps;
            temp1=lambda(i+1,1);
            temp2=kx(i+1,1);
            temp3=ky(i+1,1);
            if i+1<COL
                temp4=cnstx(i+1,1);
            end
            if i<COL
                temp5=cnstx(i,1);
            end
            temp6=krg(i+1,1);
            temp7=krw(i+1,1);
            temp8=pcgw(i+1,1);

            if lam_on==1
                lambda(i+1,1)=(1-phi(i+1,1))*lambda_s+phi(i+1,1)*(sw(i+1,1)*lambda_w+(1-sw(i+1,1)-si(i+1,1))*lambda_g+si(i+1,1)*lambda_i);
            end

            if sh(i+1,1)>0 || si(i+1,1)>0
                tmp=(1-si(i+1,1)-sh(i+1,1))^2;
            else
                tmp=1;
            end
            kx(i+1,1)=k_temp(i+1,1)*tmp;
            ky(i+1,1)=k_temp(i+1,1)*tmp;

            if (i+1) < COL
                cnstx(i+1,1)=2.0 * dz(i+2,1)*dy(i+2,1)*dz(i+1,1)*dy(i+1,1)*kx(i+2,1) * kx(i+1,1)/ ...
                    (dx(i+1,1)*dz(i+2,1)*dy(i+2,1)*kx(i+2,1) + dx(i+2,1)*dz(i+1,1)*dy(i+1,1)*kx(i+1,1));
            end
            if i < COL
                cnstx(i,1)=2.0 * dz(i+1,1)*dy(i+1,1)*dz(i,1)*dy(i,1)*kx(i+1,1) * kx(i,1)/ ...
                    (dx(i,1)*dz(i+1,1)*dy(i+1,1)*kx(i+1,1) + dx(i+1,1)*dz(i,1)*dy(i,1)*kx(i,1));
            end


            if(1-sw(i+1,1)-si(i+1,1))/(1-sh(i+1,1)-si(i+1,1)) >= sgr
                krg(i+1,1)=((1-sw(i+1,1)-si(i+1,1))/(1-sh(i+1,1)-si(i+1,1))-sgr)^gn;
            else
                krg(i+1,1)=0;
            end

            if sw(i+1,1)/(1-sh(i+1,1)-si(i+1,1))>swr
                krw(i+1,1)=(sw(i+1,1)/(1-sh(i+1,1)-si(i+1,1))-swr)^wn;
            else
                krw(i+1,1)=0;
            end

            if sh(i+1,1)>0 || si(i+1,1)>0
                tmp_ratio=sqrt(1/(1+si(i+1,1)+sh(i+1,1)+2*(1-si(i+1,1)-sh(i+1,1))/log(si(i+1,1)+sh(i+1,1))));
            else
                tmp_ratio=1;
            end
            pcgw(i+1,1)=Pd0(i+1,1)*tmp_ratio*...
                interp(tsg,tpcgw,nrpermg,(1-sw(i+1,1)-si(i+1,1))/(1-sh(i+1,1)-si(i+1,1)),1);


        else

            %*****************************************************************
            %********change si to sh as primary viarable**********************
            %*****************************************************************

            sh(i+1,1)=sh(i+1,1)+eps;
            temp1=lambda(i+1,1);
            temp2=kx(i+1,1);
            temp3=ky(i+1,1);
            if i+1<COL
                temp4=cnstx(i+1,1);
            end

            if i<COL
                temp5=cnstx(i,1);
            end
            temp6=krg(i+1,1);
            temp7=krw(i+1,1);
            temp8=pcgw(i+1,1);

            if lam_on==1
                if INDC1(i+1,1)==1
                    lambda(i+1,1)=(1-phi(i+1,1))*lambda_s+phi(i+1,1)*(sw(i+1,1)*lambda_w+sh(i+1,1)*lambda_h+(1-sh(i+1,1)-sw(i+1,1))*lambda_i);
                else
                    lambda(i+1,1)=(1-phi(i+1,1))*lambda_s+phi(i+1,1)*(sw(i+1,1)*lambda_w+(1-sw(i+1,1)-sh(i+1,1)-si(i+1,1))*lambda_g+sh(i+1,1)*lambda_h+si(i+1,1)*lambda_i);
                end
            end
            if sh(i+1,1)>0 || si(i+1,1)>0
                tmp=(1-si(i+1,1)-sh(i+1,1))^2;
            else
                tmp=1;
            end
            kx(i+1,1)=k_temp(i+1,1)*tmp;
            ky(i+1,1)=k_temp(i+1,1)*tmp;

            if (i+1) < COL
                cnstx(i+1,1)=2.0 * dz(i+2,1)*dy(i+2,1)*dz(i+1,1)*dy(i+1,1)*kx(i+2,1) * kx(i+1,1)/ ...
                    (dx(i+1,1)*dz(i+2,1)*dy(i+2,1)*kx(i+2,1) + dx(i+2,1)*dz(i+1,1)*dy(i+1,1)*kx(i+1,1));
            end
            if i < COL
                cnstx(i,1)=2.0 * dz(i+1,1)*dy(i+1,1)*dz(i,1)*dy(i,1)*kx(i+1,1) * kx(i,1)/ ...
                    (dx(i,1)*dz(i+1,1)*dy(i+1,1)*kx(i+1,1) + dx(i+1,1)*dz(i,1)*dy(i,1)*kx(i,1));
            end

            if sw(i+1,1)/(1-sh(i+1,1)-si(i+1,1))>swr
                krw(i+1,1)=(sw(i+1,1)/(1-sh(i+1,1)-si(i+1,1))-swr)^wn;
            else
                krw(i+1,1)=0;
            end

            if INDC1(i+1,1)==1
                krg(i+1,1)=0;
            else
                if (1-sw(i+1,1)-sh(i+1,1)-si(i+1,1))/(1-sh(i+1,1)-si(i+1,1)) >= sgr
                    krg(i+1,1)=((1-sw(i+1,1)-sh(i+1,1)-si(i+1,1))/(1-sh(i+1,1)-si(i+1,1))-sgr)^gn;
                else
                    krg(i+1,1)=0;
                end
            end

            if sh(i+1,1)>0 || si(i+1,1)>0
                tmp_ratio=sqrt(1/(1+si(i+1,1)+sh(i+1,1)+2*(1-si(i+1,1)-sh(i+1,1))/log(si(i+1,1)+sh(i+1,1))));
            else
                tmp_ratio=1;
            end
            pcgw(i+1,1)=Pd0(i+1,1)*tmp_ratio*...
                interp(tsg,tpcgw,nrpermg,(1-sw(i+1,1)-sh(i+1,1)-si(i+1,1))/(1-sh(i+1,1)-si(i+1,1)),1);
        end

        calc_trans(i,1);
        rsidw1= calc_Rw(i,1);
        rsidg1= calc_Rg(i,1);
        rsidc1=calc_Rc(i,1);
        rsidt1=calc_Rt(i,1);
        jac(l-3,l+4*ROW-1)=(rsidw1-rsidw0)/eps;
        jac(l-2,l+4*ROW-1)=(rsidg1-rsidg0)/eps;
        jac(l-1,l+4*ROW-1)=(rsidc1-rsidc0)/eps;
        jac(l,l+4*ROW-1)=(rsidt1-rsidt0)/eps;

        if INDC2(i+1,1) <= 3 && INDC3(i+1,1)<3
            cl(i+1,1)=cl_dimensionless(i+1,1)*cl_scale;
            dnw(i+1,1)=temp1;
        elseif INDC3(i+1,1)>2 && INDC2(i+1,1)<4

            %*****************************************************************
            %********change cl to cm as primary viarable**********************
            %*****************************************************************

            cm(i+1,1)=cm_dimensionless(i+1,1)*cm_scale;
            dnw(i+1,1)=temp1;

        elseif INDC3(i+1,1)>2 && INDC1(i+1,1)==0 && INDC2(i+1,1)<6

            %*****************************************************************
            %********change cm to si as primary viarable**********************
            %*****************************************************************

            si(i+1,1)=si(i+1,1)-eps;
            lambda(i+1,1)=temp1;
            kx(i+1,1)=temp2;
            ky(i+1,1)=temp3;
            if i+1<COL
                cnstx(i+1,1)=temp4;
            end
            if i<COL
                cnstx(i,1)=temp5;
            end
            krg(i+1,1)=temp6;
            krw(i+1,1)=temp7;
            pcgw(i+1,1)=temp8;


        else

            %*****************************************************************
            %********change si to sh as primary viarable**********************
            %*****************************************************************

            sh(i+1,1)=sh(i+1,1)-eps;
            lambda(i+1,1)=temp1;
            kx(i+1,1)=temp2;
            ky(i+1,1)=temp3;
            if i+1<COL
                cnstx(i+1,1)=temp4;
            end
            if i<COL
                cnstx(i,1)=temp5;
            end
            krg(i+1,1)=temp6;
            krw(i+1,1)=temp7;
            pcgw(i+1,1)=temp8;
        end

        calc_trans(i,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%          Change T                %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T(i+1,1)=(T_dimensionless(i+1,1)+eps_T)*T_scale+T_ini(i+1,1);
        temp1=cl(i+1,1);
        temp2=cm(i+1,1);
        temp3=dnw(i+1,1);
        temp4=dng(i+1,1);

        if INDC2(i+1,1) > 3 && INDC3(i+1,1)==1
            cl(i+1,1)=interpolation2(p_salinity,T_salinity,salinity,pw(i+1,1)/1e6,T(i+1,1));
        end
        if INDC3(i+1,1)>2
            cl(i+1,1)=((49.462^2-4*164.49*T(i+1,1))^0.5-49.462)/(2*164.49);
        end

        if (INDC2(i+1,1) > 1 && INDC3(i+1,1)==1) || (INDC3(i+1,1)>2 && INDC2(i+1,1)>3)
            if INDC1(i+1,1) == 0
                cm(i+1,1)=methane_solubility(pw(i+1,1)/1e6,T(i+1,1),2*cl(i+1,1)/0.11099/(1-cl(i+1,1)))*0.016;
            else
                cm(i+1,1)=hydrate_solubility(pw(i+1,1)/1e6,T(i+1,1),2*cl(i+1,1)/0.11099/(1-cl(i+1,1)))*0.016;
            end
        end

        dnw(i+1,1)=brine_density(pw(i+1,1)/1e6,T(i+1,1),2*cl(i+1,1)/0.11099/(1-cl(i+1,1)),...
            cm(i+1,1)/0.016);
        dng(i+1,1)=gas_density(pw(i+1,1)/1e6,T(i+1,1));

        calc_trans(i,1);
        rsidw1= calc_Rw(i,1);
        rsidg1= calc_Rg(i,1);
        rsidc1=calc_Rc(i,1);
        rsidt1=calc_Rt(i,1);
        jac(l-3,l+4*ROW)=(rsidw1-rsidw0)/eps_T;
        jac(l-2,l+4*ROW)=(rsidg1-rsidg0)/eps_T;
        jac(l-1,l+4*ROW)=(rsidc1-rsidc0)/eps_T;
        jac(l,l+4*ROW)=(rsidt1-rsidt0)/eps_T;

        T(i+1,1)=T_dimensionless(i+1,1)*T_scale+T_ini(i+1,1);
        cl(i+1,1)=temp1;
        cm(i+1,1)=temp2;
        dnw(i+1,1)=temp3;
        dng(i+1,1)=temp4;
        calc_trans(i,1);

    end   % end of if i < COL

end % end of for i=1:COL


%% ------------- Update Primary Variable ----------------
jac_temp=zeros(1,(COL-1)*ROW*4);
aa=1;
for i=2:COL
    jac_temp(aa:aa+3)=[((i-1)*ROW+j-1)*4+1 ((i-1)*ROW+j-1)*4+2 ((i-1)*ROW+j-1)*4+3  ((i-1)*ROW+j-1)*4+4];
    aa=aa+4;
end

delta=real(jac(jac_temp,jac_temp))\real(B(jac_temp));

l=0;
for i=2:COL
    for j=1:ROW
        l=l+4;
        if INDC2(i,1) == 1
            dpw=delta(l-3);
            dcm=delta(l-2);
            dcl=delta(l-1);
            dT=delta(l);

            pw(i,1)=(pw_dimensionless(i,1)+dpw)*pw_scale+pw_ini(i,1);
            cm(i,1)=(cm_dimensionless(i,1)+dcm)*cm_scale;
            cl(i,1)=(cl_dimensionless(i,1)+dcl)*cl_scale;
            T(i,1)=(T_dimensionless(i,1)+dT)*T_scale+T_ini(i,1);

            sw(i,1)=1.0;
            sg(i,1)=0;
            sh(i,1)=0;
            si(i,1)=0;
        elseif INDC2(i,1) == 2 || INDC2(i,1) == 3
            if INDC3(i,1)==1
                dpw=delta(l-3);
                dsw=delta(l-2);
                dcl=delta(l-1);
                dT=delta(l);

                pw(i,1)=(pw_dimensionless(i,1)+dpw)*pw_scale+pw_ini(i,1);
                sw(i,1)=sw(i,1)+dsw;
                cl(i,1)=(cl_dimensionless(i,1)+dcl)*cl_scale;
                T(i,1)=(T_dimensionless(i,1)+dT)*T_scale+T_ini(i,1);

                if INDC1(i,1)==0
                    sg(i,1)=1-sw(i,1);
                    sh(i,1)=0;
                    si(i,1)=0;
                else
                    sh(i,1)=1-sw(i,1);
                    sg(i,1)=0;
                    si(i,1)=0;
                end

            elseif INDC3(i,1)==2
                dpw=delta(l-3);
                dsi=delta(l-2);
                dcl=delta(l-1);
                dT=delta(l);

                pw(i,1)=(pw_dimensionless(i,1)+dpw)*pw_scale+pw_ini(i,1);
                si(i,1)=si(i,1)+dsi;
                cl(i,1)=(cl_dimensionless(i,1)+dcl)*cl_scale;
                T(i,1)=(T_dimensionless(i,1)+dT)*T_scale+T_ini(i,1);

                if INDC1(i,1)==0
                    sg(i,1)=1-si(i,1);
                    sh(i,1)=0;
                    sw(i,1)=0;
                    cm(i,1)=0;
                    cl(i,1)=0;
                else
                    sh(i,1)=1-si(i,1);
                    sg(i,1)=0;
                    sw(i,1)=0;
                    cm(i,1)=0;
                    cl(i,1)=0;
                end

            else

                dpw=delta(l-3);
                dsw=delta(l-2);
                dcm=delta(l-1);
                dT=delta(l);

                pw(i,1)=(pw_dimensionless(i,1)+dpw)*pw_scale+pw_ini(i,1);
                sw(i,1)=sw(i,1)+dsw;
                cm(i,1)=(cm_dimensionless(i,1)+dcm)*cm_scale;
                T(i,1)=(T_dimensionless(i,1)+dT)*T_scale+T_ini(i,1);


                si(i,1)=1-sw(i,1);
                sh(i,1)=0;
                sg(i,1)=0;


            end

        elseif INDC2(i,1) == 4 || INDC2(i,1) == 5
            if INDC3(i,1)==1
                dpw=delta(l-3);
                dsw=delta(l-2);
                dsh=delta(l-1);
                dT=delta(l);

                pw(i,1)=(pw_dimensionless(i,1)+dpw)*pw_scale+pw_ini(i,1);
                sw(i,1)=sw(i,1)+dsw;
                sh(i,1)=sh(i,1)+dsh;
                T(i,1)=(T_dimensionless(i,1)+dT)*T_scale+T_ini(i,1);


                sg(i,1)=1-sw(i,1)-sh(i,1);
                si(i,1)=0;

            elseif INDC3(i,1)==2
                dpw=delta(l-3);
                dsi=delta(l-2);
                dsh=delta(l-1);
                dT=delta(l);

                pw(i,1)=(pw_dimensionless(i,1)+dpw)*pw_scale+pw_ini(i,1);
                si(i,1)=si(i,1)+dsi;
                sh(i,1)=sh(i,1)+dsh;
                T(i,1)=(T_dimensionless(i,1)+dT)*T_scale+T_ini(i,1);


                sg(i,1)=1-si(i,1)-sh(i,1);
                sw(i,1)=0;
                cl(i,1)=0;
                cm(i,1)=0;

            elseif INDC3(i,1)>2 && INDC1(i,1)==1
                dpw=delta(l-3);
                dsw=delta(l-2);
                dsh=delta(l-1);
                dT=delta(l);

                pw(i,1)=(pw_dimensionless(i,1)+dpw)*pw_scale+pw_ini(i,1);
                sw(i,1)=sw(i,1)+dsw;
                sh(i,1)=sh(i,1)+dsh;
                T(i,1)=(T_dimensionless(i,1)+dT)*T_scale+T_ini(i,1);


                si(i,1)=1-sw(i,1)-sh(i,1);
                sg(i,1)=0;

            else
                dpw=delta(l-3);
                dsw=delta(l-2);
                dsi=delta(l-1);
                dT=delta(l);

                pw(i,1)=(pw_dimensionless(i,1)+dpw)*pw_scale+pw_ini(i,1);
                sw(i,1)=sw(i,1)+dsw;
                si(i,1)=si(i,1)+dsi;
                T(i,1)=(T_dimensionless(i,1)+dT)*T_scale+T_ini(i,1);


                sg(i,1)=1-sw(i,1)-si(i,1);
                sh(i,1)=0;

            end

        else
            dsi=delta(l-3);
            dsw=delta(l-2);
            dsh=delta(l-1);
            dT=delta(l);

            si(i,1)=si(i,1)+dsi;
            sw(i,1)=sw(i,1)+dsw;
            sh(i,1)=sh(i,1)+dsh;
            T(i,1)=(T_dimensionless(i,1)+dT)*T_scale+T_ini(i,1);

            sg(i,1)=1-sw(i,1)-si(i,1)-sh(i,1);
        end % end of if INDC2(i,1)==1
    end % end of for j=1:ROW
end % end of for i=1:COL


