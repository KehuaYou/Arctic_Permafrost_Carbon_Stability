global COL ROW
global INDC2 INDC3 INDC1
global dx dy dz
global pw  T 
global cl cm
global sw sg sh si 
global salinity
global dnw dng
global k_temp kx  cnstx 
global krw krg sgr swr wn gn 
global lambda lambda_g lambda_h lambda_w lambda_s lambda_i Pd0 pcgw phi
global tsg tpcgw nrpermg


for i=2:COL
    for j=1:ROW
        
        if INDC2(i,j) > 1
            %-------- update salinity ---------------
            if INDC2(i,j) > 3 && INDC3(i,j)==1
                cl(i,j)=interpolation2(p_salinity,T_salinity,salinity,pw(i,j)/1e6,T(i,j));
            elseif INDC3(i,j)>2
                cl(i,j)=((49.462^2-4*164.49*T(i,j))^0.5-49.462)/(2*164.49);
            elseif  INDC3(i,1)==2
                cl(i,j)=0;
            end
            %-------- update methane solubility ---------------
            if INDC3(i,j)==2
                cm(i,j)=0;
            elseif INDC3(i,1)==1 || (INDC3(i,1)>2 && INDC2(i,1)>3)
                if INDC1(i,j) == 0
                    cm(i,j)=methane_solubility(pw(i,j)/1e6,T(i,j),2*cl(i,j)/0.11099/(1-cl(i,j)))*0.016;
                else
                    cm(i,j)=hydrate_solubility(pw(i,j)/1e6,T(i,j),2*cl(i,j)/0.11099/(1-cl(i,j)))*0.016;
                end
            end
            %-------- update salinity, pressure and methane solubility ---------------
            if INDC2(i,j)>=6
                    cl(i,j)=((49.462^2-4*164.49*T(i,j))^0.5-49.462)/(2*164.49);
                    pw(i,j)=interpolation_pressure(p_salinity,T_salinity,salinity,cl(i,j),T(i,j))*1e6;
                    cm(i,j)=methane_solubility(pw(i,j)/1e6,T(i,j),2*cl(i,j)/0.11099/(1-cl(i,j)))*0.016;
            end 
        end

        %-------- update densities ---------------
        dnw(i,j)=brine_density(pw(i,j)/1e6,T(i,j),2*cl(i,j)/0.11099/(1-cl(i,j)),cm(i,j)/0.016);
        dng(i,j)=gas_density(pw(i,j)/1e6,T(i,j));

        %-------- update bulk thermal conductivity ---------------
        lambda(i,j)=(1-phi(i,j))*lambda_s+phi(i,j)*(sw(i,j)*lambda_w+sg(i,j)*lambda_g+sh(i,j)*lambda_h+si(i,j)*lambda_i);

        %-------- update sediment effective permeability ---------------
        if sh(i,1)>0 || si(i,1)>0
            tmp=(1-si(i,1)-sh(i,1))^2;
        else
            tmp=1;
        end
        kx(i,1)=k_temp(i,1)*tmp;
        %-------- update relative permeability ---------------
        if sw(i,j)/(1-sh(i,j)-si(i,j))>swr
            krw(i,j)=(sw(i,j)/(1-sh(i,j)-si(i,j))-swr)^wn;
        else
            krw(i,j)=0;
        end
        if sg(i,j)/(1-sh(i,j)-si(i,j)) > sgr
            krg(i,j)=(sg(i,j)/(1-sh(i,j)-si(i,j))-sgr)^gn;
        else
            krg(i,j)=0;
        end
        
        %-------- update gas-water capillary pressure ---------------            
        if sh(i,j)>0 || si(i,j)>0
            tmp_ratio=sqrt(1/(1+si(i,1)+sh(i,1)+2*(1-si(i,1)-sh(i,1))/log(si(i,1)+sh(i,1))));
        else
            tmp_ratio=1;
        end
        pcgw(i,j)=Pd0(i,j)*tmp_ratio*interp(tsg,tpcgw,nrpermg,sg(i,j)/(1-sh(i,j)-si(i,j)),1);
    end
end
for i=1:COL
    for j=1:ROW
        if i < COL
            cnstx(i,j)=2.0 * dz(i+1,j)*dy(i+1,j)*dz(i,j)*dy(i,j)*kx(i+1,j) * kx(i,j)/ ...
                (dx(i,j)*dz(i+1,j)*dy(i+1,j)*kx(i+1,j) + dx(i+1,j)*dz(i,j)*dy(i,j)*kx(i,j));
        else
            cnstx(i,j)=0.0;
        end
    end
end

