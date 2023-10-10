function Rg=calc_Rg(i,j)
global COL
global pw  pcgw  dpth
global tgx tgx1  twx twx1
global dnw dnw_0 sw sw_0 dng dng_0 sg_0 sh sh_0 cl cl_0 cm cm_0 si
global phi phi_0 dt
global dnh Mch4
global INDC1 INDC2 INDC3 
global vb
global qg
global dy dz
global qg_biogenic
global pw_top cm_top cl_top
global Dmethanex


Rg=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if i==1
    
    if twx(i,j)*(pw(i+1,j)-pw(i,j)) > twx1(i,j)*(dpth(i+1,j)-dpth(i,j))
        Rg=Rg + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i+1,j)*(1-cl(i+1,j))*cm(i+1,j);
        Rg=Rg - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i+1,j)*(1-cl(i+1,j))*cm(i+1,j);
    else
        Rg=Rg + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i,j)*(1-cl(i,j))*cm(i,j);
        Rg=Rg - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i,j)*(1-cl(i,j))*cm(i,j);
    end
    if tgx(i,j)*(pw(i+1,j)+pcgw(i+1,j)-pw(i,j)-pcgw(i,j)) > tgx1(i,j)*(dpth(i+1,j)-dpth(i,j))
        Rg=Rg + tgx(i,j)*(pw(i+1,j)-pw(i,j))*dng(i+1,j) + tgx(i,j)*(pcgw(i+1,j)-pcgw(i,j))*dng(i+1,j);
        Rg=Rg - tgx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dng(i+1,j);
    else
        Rg=Rg + tgx(i,j)*(pw(i+1,j)-pw(i,j))*dng(i,j) + tgx(i,j)*(pcgw(i+1,j)-pcgw(i,j))*dng(i,j);
        Rg=Rg - tgx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dng(i,j);
    end
    Rg=Rg + Dmethanex(i,j)*(phi(i+1,j)*sw(i+1,j)+phi(i,j)*sw(i,j))/2*(dnw(i+1,j)+dnw(i,j))/2*(cm(i+1,j)-cm(i,j));
    %++++++++++++++++++++++++++fixed gas pressure++++++++++++++++++++++++++++++
    
    if twx(i,j)*(pw(i,j)-pw_top) > twx1(i,j)*(dpth(i,j)-0)
        Rg=Rg - 2*twx(i,j)*(pw(i,j)-pw_top)*dnw(i,j)*(1-cl(i,j))*cm(i,j);
        Rg=Rg + 2*twx1(i,j)*(dpth(i,j)-0)*dnw(i,j)*(1-cl(i,j))*cm(i,j);
    else
        Rg=Rg - 2*twx(i,j)*(pw(i,j)-pw_top)*dnw(i,j)*(1-cl_top)*cm_top;
        Rg=Rg + 2*twx1(i,j)*(dpth(i,j)-0)*dnw(i,j)*(1-cl_top)*cm_top;
    end
    
    if tgx(i,j)*(pw(i,j)+pcgw(i,j)-pw_top) > tgx1(i,j)*(dpth(i,j)-0)
        Rg=Rg - 2*tgx(i,j)*(pw(i,j)+pcgw(i,j)-pw_top)*dng(i,j);
        Rg=Rg + 2*tgx1(i,j)*(dpth(i,j)-0)*dng(i,j);
    else
        Rg=Rg - 2*tgx(i,j)*(pw(i,j)+pcgw(i,j)-pw_top)*dng(i,j);
        Rg=Rg + 2*tgx1(i,j)*(dpth(i,j)-0)*dng(i,j);
    end
    Rg=Rg - 2*Dmethanex(i,j)*phi(i,j)*sw(i,j)*dnw(i,j)*(cm(i,j)-cm_top);
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
elseif i==COL
    
    if twx(i-1,j)*(pw(i,j)-pw(i-1,j)) > twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
        Rg=Rg - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i,j)*(1-cl(i,j))*cm(i,j);
        Rg=Rg + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i,j)*(1-cl(i,j))*cm(i,j);
    else
        Rg=Rg - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i-1,j)*(1-cl(i-1,j))*cm(i-1,j);
        Rg=Rg + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i-1,j)*(1-cl(i-1,j))*cm(i-1,j);
    end
    if tgx(i-1,j)*(pw(i,j)+pcgw(i,j)-pw(i-1,j)-pcgw(i-1,j)) > tgx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
        Rg=Rg - tgx(i-1,j)*(pw(i,j)-pw(i-1,j))*dng(i,j) - tgx(i-1,j)*(pcgw(i,j)-pcgw(i-1,j))*dng(i,j);
        Rg=Rg + tgx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dng(i,j);
    else
        Rg=Rg - tgx(i-1,j)*(pw(i,j)-pw(i-1,j))*dng(i-1,j) - tgx(i-1,j)*(pcgw(i,j)-pcgw(i-1,j))*dng(i-1,j);
        Rg=Rg + tgx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dng(i-1,j);
    end
    Rg=Rg - Dmethanex(i-1,j)*(phi(i,j)*sw(i,j)+phi(i-1,j)*sw(i-1,j))/2*(dnw(i,j)+dnw(i-1,j))/2*(cm(i,j)-cm(i-1,j));
    
else
    if twx(i,j)*(pw(i+1,j)-pw(i,j)) > twx1(i,j)*(dpth(i+1,j)-dpth(i,j))
        Rg=Rg + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i+1,j)*(1-cl(i+1,j))*cm(i+1,j);
        Rg=Rg - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i+1,j)*(1-cl(i+1,j))*cm(i+1,j);
    else
        Rg=Rg + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i,j)*(1-cl(i,j))*cm(i,j);
        Rg=Rg - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i,j)*(1-cl(i,j))*cm(i,j);
    end
    if tgx(i,j)*(pw(i+1,j)+pcgw(i+1,j)-pw(i,j)-pcgw(i,j)) > tgx1(i,j)*(dpth(i+1,j)-dpth(i,j))
        Rg=Rg + tgx(i,j)*(pw(i+1,j)-pw(i,j))*dng(i+1,j) + tgx(i,j)*(pcgw(i+1,j)-pcgw(i,j))*dng(i+1,j);
        Rg=Rg - tgx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dng(i+1,j);
    else
        Rg=Rg + tgx(i,j)*(pw(i+1,j)-pw(i,j))*dng(i,j) + tgx(i,j)*(pcgw(i+1,j)-pcgw(i,j))*dng(i,j);
        Rg=Rg - tgx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dng(i,j);
    end
    Rg=Rg + Dmethanex(i,j)*(phi(i+1,j)*sw(i+1,j)+phi(i,j)*sw(i,j))/2*(dnw(i+1,j)+dnw(i,j))/2*(cm(i+1,j)-cm(i,j));
    
    if twx(i-1,j)*(pw(i,j)-pw(i-1,j)) > twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
        Rg=Rg - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i,j)*(1-cl(i,j))*cm(i,j);
        Rg=Rg + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i,j)*(1-cl(i,j))*cm(i,j);
    else
        Rg=Rg - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i-1,j)*(1-cl(i-1,j))*cm(i-1,j);
        Rg=Rg + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i-1,j)*(1-cl(i-1,j))*cm(i-1,j);
    end
    if tgx(i-1,j)*(pw(i,j)+pcgw(i,j)-pw(i-1,j)-pcgw(i-1,j)) > tgx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
        Rg=Rg - tgx(i-1,j)*(pw(i,j)-pw(i-1,j))*dng(i,j) - tgx(i-1,j)*(pcgw(i,j)-pcgw(i-1,j))*dng(i,j);
        Rg=Rg + tgx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dng(i,j);
    else
        Rg=Rg - tgx(i-1,j)*(pw(i,j)-pw(i-1,j))*dng(i-1,j) - tgx(i-1,j)*(pcgw(i,j)-pcgw(i-1,j))*dng(i-1,j);
        Rg=Rg + tgx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dng(i-1,j);
    end
    
    Rg=Rg - Dmethanex(i-1,j)*(phi(i,j)*sw(i,j)+phi(i-1,j)*sw(i-1,j))/2*(dnw(i,j)+dnw(i-1,j))/2*(cm(i,j)-cm(i-1,j));
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if INDC3(i,j)==1
    if INDC2(i,j) == 1
        
        Rg=-Rg+vb(i,j)*((phi(i,j)*1*dnw(i,j)*(1-cl(i,j))*cm(i,j)+phi(i,j)*0*dnh*Mch4+phi(i,j)*0*dng(i,j))-...
            (phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)+phi_0(i,j)*sh_0(i,j)*dnh*Mch4+phi_0(i,j)*sg_0(i,j)*dng_0(i,j)))/dt;
        
    elseif INDC2(i,j) == 2
        
        if INDC1(i,j) == 0
            
            Rg=-Rg+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))*cm(i,j)+phi(i,j)*0*dnh*Mch4+phi(i,j)*(1-sw(i,j))*dng(i,j))-...
            (phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)+phi_0(i,j)*sh_0(i,j)*dnh*Mch4+phi_0(i,j)*sg_0(i,j)*dng_0(i,j)))/dt;
        
        else
            
            Rg=-Rg+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))*cm(i,j)+phi(i,j)*(1-sw(i,j))*dnh*Mch4+phi(i,j)*0*dng(i,j))-...
            (phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)+phi_0(i,j)*sh_0(i,j)*dnh*Mch4+phi_0(i,j)*sg_0(i,j)*dng_0(i,j)))/dt;
            
        end
        
    elseif INDC2(i,j) == 3
        
        if INDC1(i,j) == 0
            
            Rg=-Rg+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))*cm(i,j)+phi(i,j)*0*dnh*Mch4+phi(i,j)*(1-sw(i,j))*dng(i,j))-...
            (phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)+phi_0(i,j)*sh_0(i,j)*dnh*Mch4+phi_0(i,j)*sg_0(i,j)*dng_0(i,j)))/dt;
        
          
        else
            
            Rg=-Rg+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))*cm(i,j)+phi(i,j)*(1-sw(i,j))*dnh*Mch4+phi(i,j)*0*dng(i,j))-...
            (phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)+phi_0(i,j)*sh_0(i,j)*dnh*Mch4+phi_0(i,j)*sg_0(i,j)*dng_0(i,j)))/dt;
            
        end
    elseif INDC2(i,j) == 4 || INDC2(i,j) == 5
        
        Rg=-Rg+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))*cm(i,j)+phi(i,j)*sh(i,j)*dnh*Mch4+phi(i,j)*(1-sw(i,j)-sh(i,j))*dng(i,j))-...
            (phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)+phi_0(i,j)*sh_0(i,j)*dnh*Mch4+phi_0(i,j)*sg_0(i,j)*dng_0(i,j)))/dt;
            
    else
        
        Rg=-Rg+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))*cm(i,j)+phi(i,j)*sh(i,j)*dnh*Mch4+phi(i,j)*(1-sw(i,j)-sh(i,j)-si(i,j))*dng(i,j))-...
            (phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)+phi_0(i,j)*sh_0(i,j)*dnh*Mch4+phi_0(i,j)*sg_0(i,j)*dng_0(i,j)))/dt;

    end
end




if INDC3(i,j)==2
    if INDC2(i,j) == 2 || INDC2(i,j) ==3
        if INDC1(i,j) == 0
            
            Rg=-Rg+vb(i,j)*((phi(i,j)*0*dnw(i,j)*(1-cl(i,j))*cm(i,j)+phi(i,j)*0*dnh*Mch4+phi(i,j)*(1-si(i,j))*dng(i,j))-...
            (phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)+phi_0(i,j)*sh_0(i,j)*dnh*Mch4+phi_0(i,j)*sg_0(i,j)*dng_0(i,j)))/dt;

        else
           
            Rg=-Rg+vb(i,j)*((phi(i,j)*0*dnw(i,j)*(1-cl(i,j))*cm(i,j)+phi(i,j)*(1-si(i,j))*dnh*Mch4+phi(i,j)*0*dng(i,j))-...
            (phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)+phi_0(i,j)*sh_0(i,j)*dnh*Mch4+phi_0(i,j)*sg_0(i,j)*dng_0(i,j)))/dt;

        end
        
    else
        
        Rg=-Rg+vb(i,j)*((phi(i,j)*0*dnw(i,j)*(1-cl(i,j))*cm(i,j)+phi(i,j)*sh(i,j)*dnh*Mch4+phi(i,j)*(1-si(i,j)-sh(i,j))*dng(i,j))-...
            (phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)+phi_0(i,j)*sh_0(i,j)*dnh*Mch4+phi_0(i,j)*sg_0(i,j)*dng_0(i,j)))/dt;

    end
end





if INDC3(i,j)==3
    if INDC2(i,j)==2 % L-> L+I
        
         Rg=-Rg+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))*cm(i,j)+phi(i,j)*0*dnh*Mch4+phi(i,j)*0*dng(i,j))-...
            (phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)+phi_0(i,j)*sh_0(i,j)*dnh*Mch4+phi_0(i,j)*sg_0(i,j)*dng_0(i,j)))/dt;

    elseif INDC2(i,j)==4
     
            if INDC1(i,j)==0    % L+G -> L+G+I, I+G -> I+G+L
                
                Rg=-Rg+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))*cm(i,j)+phi(i,j)*0*dnh*Mch4+phi(i,j)*(1-sw(i,j)-si(i,j))*dng(i,j))-...
            (phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)+phi_0(i,j)*sh_0(i,j)*dnh*Mch4+phi_0(i,j)*sg_0(i,j)*dng_0(i,j)))/dt;

            else    % L+H -> L+H+I, I+H -> I+H+L
                
                Rg=-Rg+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))*cm(i,j)+phi(i,j)*sh(i,j)*dnh*Mch4+phi(i,j)*0*dng(i,j))-...
            (phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)+phi_0(i,j)*sh_0(i,j)*dnh*Mch4+phi_0(i,j)*sg_0(i,j)*dng_0(i,j)))/dt;

            end
        
        
    elseif INDC2(i,j)==6  
        
         Rg=-Rg+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))*cm(i,j)+phi(i,j)*sh(i,j)*dnh*Mch4+phi(i,j)*(1-sw(i,j)-sh(i,j)-si(i,j))*dng(i,j))-...
            (phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)+phi_0(i,j)*sh_0(i,j)*dnh*Mch4+phi_0(i,j)*sg_0(i,j)*dng_0(i,j)))/dt;

    end
end






if INDC3(i,j)==4
    if INDC2(i,j)==3   % L+I
        
        Rg=-Rg+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))*cm(i,j)+phi(i,j)*0*dnh*Mch4+phi(i,j)*0*dng(i,j))-...
            (phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)+phi_0(i,j)*sh_0(i,j)*dnh*Mch4+phi_0(i,j)*sg_0(i,j)*dng_0(i,j)))/dt;

    elseif INDC2(i,j)==4 || INDC2(i,j)==5
        
        if INDC1(i,j)==0 % L+I -> L+I+G, L+I+H+G -> L+I+G
           Rg=-Rg+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))*cm(i,j)+phi(i,j)*0*dnh*Mch4+phi(i,j)*(1-sw(i,j)-si(i,j))*dng(i,j))-...
            (phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)+phi_0(i,j)*sh_0(i,j)*dnh*Mch4+phi_0(i,j)*sg_0(i,j)*dng_0(i,j)))/dt;
 
        else   % L+I -> L+I+H, L+I+H+G-> L+I+H
            Rg=-Rg+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))*cm(i,j)+phi(i,j)*sh(i,j)*dnh*Mch4+phi(i,j)*0*dng(i,j))-...
            (phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)+phi_0(i,j)*sh_0(i,j)*dnh*Mch4+phi_0(i,j)*sg_0(i,j)*dng_0(i,j)))/dt;
            
        end
    elseif INDC2(i,j)==6 || INDC2(i,j)==7  % L+I+H-> L+I+H+G, L+I+G-> L+I+H+G
        
        Rg=-Rg+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))*cm(i,j)+phi(i,j)*sh(i,j)*dnh*Mch4+phi(i,j)*(1-sw(i,j)-sh(i,j)-si(i,j))*dng(i,j))-...
            (phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))*cm_0(i,j)+phi_0(i,j)*sh_0(i,j)*dnh*Mch4+phi_0(i,j)*sg_0(i,j)*dng_0(i,j)))/dt;
            
    end
end

Rg=Rg-vb(i,j)*qg_biogenic(i,j); % biogenic methane


%++++++++++++++++++++inject gas from below+++++++++++++++++++++++++
if i==COL
    Rg=Rg-dy(i,j)*dz(i,j)*qg/(365*24*3600);
end





















