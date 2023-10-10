function Rw=calc_Rw(i,j)
global COL
global pw  dpth
global twx twx1
global sw sw_0 dnw dnw_0 cl cl_0
global phi phi_0 sh sh_0 si si_0
global dnh Mh2o dt dni
global INDC1 INDC2 INDC3 
global vb
global qw dy dz
global pw_top  cl_top

Rw=0;
if i==1
    %++++++++++++++++++original no flow+++++++++++++++++++++++++++++++++++
    if twx(i,j)*(pw(i+1,j)-pw(i,j)) > twx1(i,j)*(dpth(i+1,j)-dpth(i,j))
        Rw=Rw + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i+1,j)*(1-cl(i+1,j));
        Rw=Rw - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i+1,j)*(1-cl(i+1,j));
    else
        Rw=Rw + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i,j)*(1-cl(i,j));
        Rw=Rw - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i,j)*(1-cl(i,j));
    end
    
    %++++++++++++++++++fixed outside water pressure++++++++++++++++++++++++
    if twx(i,j)*(pw(i,j)-pw_top) > twx1(i,j)*(dpth(i,j)-0)
        Rw=Rw - 2*twx(i,j)*(pw(i,j)-pw_top)*dnw(i,j)*(1-cl(i,j));
        Rw=Rw + 2*twx1(i,j)*(dpth(i,j)-0)*dnw(i,j)*(1-cl(i,j));
    else
        Rw=Rw - 2*twx(i,j)*(pw(i,j)-pw_top)*dnw(i,j)*(1-cl_top);
        Rw=Rw + 2*twx1(i,j)*(dpth(i,j)-0)*dnw(i,j)*(1-cl_top);
    end
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
elseif i==COL
    
    if twx(i-1,j)*(pw(i,j)-pw(i-1,j)) > twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
        Rw=Rw - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i,j)*(1-cl(i,j));
        Rw=Rw + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i,j)*(1-cl(i,j));
    else
        Rw=Rw - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i-1,j)*(1-cl(i-1,j));
        Rw=Rw + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i-1,j)*(1-cl(i-1,j));
    end
   
else
    if twx(i,j)*(pw(i+1,j)-pw(i,j)) > twx1(i,j)*(dpth(i+1,j)-dpth(i,j))
        Rw=Rw + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i+1,j)*(1-cl(i+1,j));
        Rw=Rw - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i+1,j)*(1-cl(i+1,j));
    else
        Rw=Rw + twx(i,j)*(pw(i+1,j)-pw(i,j))*dnw(i,j)*(1-cl(i,j));
        Rw=Rw - twx1(i,j)*(dpth(i+1,j)-dpth(i,j))*dnw(i,j)*(1-cl(i,j));
    end
    
    if twx(i-1,j)*(pw(i,j)-pw(i-1,j)) > twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))
        Rw=Rw - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i,j)*(1-cl(i,j));
        Rw=Rw + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i,j)*(1-cl(i,j));
    else
        Rw=Rw - twx(i-1,j)*(pw(i,j)-pw(i-1,j))*dnw(i-1,j)*(1-cl(i-1,j));
        Rw=Rw + twx1(i-1,j)*(dpth(i,j)-dpth(i-1,j))*dnw(i-1,j)*(1-cl(i-1,j));
    end
end




if INDC3(i,j)==1
    if INDC2(i,j) == 1
        
        Rw=-Rw+vb(i,j)*((phi(i,j)*1*dnw(i,j)*(1-cl(i,j))+phi(i,j)*0*dni+phi(i,j)*0*dnh*Mh2o)-...
            ((phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*si_0(i,j)*dni+phi_0(i,j)*sh_0(i,j)*dnh*Mh2o)))/dt;
        
    elseif INDC2(i,j) == 2
        if INDC1(i,j) == 0
            
            Rw=-Rw+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))+phi(i,j)*0*dni+phi(i,j)*0*dnh*Mh2o)-...
            ((phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*si_0(i,j)*dni+phi_0(i,j)*sh_0(i,j)*dnh*Mh2o)))/dt;
        
        else
            
            Rw=-Rw+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))+phi(i,j)*0*dni+phi(i,j)*(1-sw(i,j))*dnh*Mh2o)-...
            ((phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*si_0(i,j)*dni+phi_0(i,j)*sh_0(i,j)*dnh*Mh2o)))/dt;
 
        end
        
    elseif INDC2(i,j) == 3
        if INDC1(i,j) == 0
            Rw=-Rw+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))+phi(i,j)*0*dni+phi(i,j)*0*dnh*Mh2o)-...
            ((phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*0*dni+phi_0(i,j)*0*dnh*Mh2o)))/dt;
 
        else
             Rw=-Rw+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))+phi(i,j)*0*dni+phi(i,j)*(1-sw(i,j))*dnh*Mh2o)-...
            ((phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*0*dni+phi_0(i,j)*sh_0(i,j)*dnh*Mh2o)))/dt;
        end
        
    elseif INDC2(i,j)==4
        
        Rw=-Rw+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))+phi(i,j)*0*dni+phi(i,j)*sh(i,j)*dnh*Mh2o)-...
            ((phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*si_0(i,j)*dni+phi_0(i,j)*sh_0(i,j)*dnh*Mh2o)))/dt;
    else
        Rw=-Rw+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))+phi(i,j)*0*dni+phi(i,j)*sh(i,j)*dnh*Mh2o)-...
            ((phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*0*dni+phi_0(i,j)*sh_0(i,j)*dnh*Mh2o)))/dt;
    end
end



if INDC3(i,j)==2
    
    if INDC2(i,j)==2
        if INDC1(i,j) == 0  % I+G+H -> I+G
            
            Rw=-Rw+vb(i,j)*((phi(i,j)*0*dnw(i,j)*(1-cl(i,j))+phi(i,j)*si(i,j)*dni+phi(i,j)*0*dnh*Mh2o)-...
            ((phi_0(i,j)*0*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*si_0(i,j)*dni+phi_0(i,j)*sh_0(i,j)*dnh*Mh2o)))/dt;
        
        else    % I+G+H -> I+H
            
            Rw=-Rw+vb(i,j)*((phi(i,j)*0*dnw(i,j)*(1-cl(i,j))+phi(i,j)*si(i,j)*dni+phi(i,j)*(1-si(i,j))*dnh*Mh2o)-...
            ((phi_0(i,j)*0*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*si_0(i,j)*dni+phi_0(i,j)*sh_0(i,j)*dnh*Mh2o)))/dt;
        
        end
    elseif INDC2(i,j) == 3
        if INDC1(i,j) == 0  % I+G 
            
            Rw=-Rw+vb(i,j)*((phi(i,j)*0*dnw(i,j)*(1-cl(i,j))+phi(i,j)*si(i,j)*dni+phi(i,j)*0*dnh*Mh2o)-...
            ((phi_0(i,j)*0*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*si_0(i,j)*dni+phi_0(i,j)*0*dnh*Mh2o)))/dt;
        
        else             % I+H
            
            Rw=-Rw+vb(i,j)*((phi(i,j)*0*dnw(i,j)*(1-cl(i,j))+phi(i,j)*si(i,j)*dni+phi(i,j)*(1-si(i,j))*dnh*Mh2o)-...
            ((phi_0(i,j)*0*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*si_0(i,j)*dni+phi_0(i,j)*sh_0(i,j)*dnh*Mh2o)))/dt;
        
        end
    elseif INDC2(i,j)==4
        if INDC1(i,j)==0   % I+G -> I+G+H
            
            Rw=-Rw+vb(i,j)*((phi(i,j)*0*dnw(i,j)*(1-cl(i,j))+phi(i,j)*si(i,j)*dni+phi(i,j)*sh(i,j)*dnh*Mh2o)-...
            ((phi_0(i,j)*0*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*si_0(i,j)*dni+phi_0(i,j)*0*dnh*Mh2o)))/dt;
        
        else               % I+H -> I+G+H
           
            Rw=-Rw+vb(i,j)*((phi(i,j)*0*dnw(i,j)*(1-cl(i,j))+phi(i,j)*si(i,j)*dni+phi(i,j)*sh(i,j)*dnh*Mh2o)-...
            ((phi_0(i,j)*0*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*si_0(i,j)*dni+phi_0(i,j)*sh_0(i,j)*dnh*Mh2o)))/dt;
        end
    else                 % I+G+H
             
        Rw=-Rw+vb(i,j)*((phi(i,j)*0*dnw(i,j)*(1-cl(i,j))+phi(i,j)*si(i,j)*dni+phi(i,j)*sh(i,j)*dnh*Mh2o)-...
            ((phi_0(i,j)*0*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*si_0(i,j)*dni+phi_0(i,j)*sh_0(i,j)*dnh*Mh2o)))/dt;
    end
end



if INDC3(i,j)==3  % INDC3(i,j)==3 only happens when Ice or Liquid water appears, only happens when phase number increases
    if INDC2(i,j) == 2 % L -> L+I
        
        Rw=-Rw+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))+phi(i,j)*(1-sw(i,j))*dni+phi(i,j)*0*dnh*Mh2o)-...
            ((phi_0(i,j)*1*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*0*dni+phi_0(i,j)*0*dnh*Mh2o)))/dt;
        
    elseif INDC2(i,j)==4
    
            if INDC1(i,j)==0    % L+G -> L+G+I, I+G -> I+G+L
                
                Rw=-Rw+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))+phi(i,j)*si(i,j)*dni+phi(i,j)*0*dnh*Mh2o)-...
            ((phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*0*dni+phi_0(i,j)*0*dnh*Mh2o)))/dt;
        
            else               % L+H -> L+H+I, I+H -> I+H+L
               
                 Rw=-Rw+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))+phi(i,j)*(1-sw(i,j)-sh(i,j))*dni+phi(i,j)*sh(i,j)*dnh*Mh2o)-...
            ((phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*0*dni+phi_0(i,j)*sh_0(i,j)*dnh*Mh2o)))/dt;
        
            end
        
        
    elseif INDC2(i,j)==6   % L+G+H -> L+G+H+I, I+G+H -> I+G+H+L
            
             Rw=-Rw+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))+phi(i,j)*si(i,j)*dni+phi(i,j)*sh(i,j)*dnh*Mh2o)-...
            ((phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*0*dni+phi_0(i,j)*sh_0(i,j)*dnh*Mh2o)))/dt;
        
    end
end




if INDC3(i,j)==4
    
    if INDC2(i,j) == 3  % L+I
        
         Rw=-Rw+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))+phi(i,j)*( 1-sw(i,j))*dni+phi(i,j)*0*dnh*Mh2o)-...
            ((phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*si_0(i,j)*dni+phi_0(i,j)*0*dnh*Mh2o)))/dt;
        
    elseif INDC2(i,j)==4 || INDC2(i,j)==5 % L+I -> L+I+H, L+I+G, or L+I+H+G -> L+I+H, L+I+G
        
        if INDC1(i,j)==0 % L+I -> L+I+G, L+I+G+H -> L+I+G
            
            Rw=-Rw+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))+phi(i,j)*si(i,j)*dni+phi(i,j)*0*dnh*Mh2o)-...
            ((phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*si_0(i,j)*dni+phi_0(i,j)*sh_0(i,j)*dnh*Mh2o)))/dt;
        
        else  % L+I -> L+I+H, L+I+G+H -> L+I+H
            
            Rw=-Rw+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))+phi(i,j)*(1-sw(i,j)-sh(i,j))*dni+phi(i,j)*sh(i,j)*dnh*Mh2o)-...
            ((phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*si_0(i,j)*dni+phi_0(i,j)*sh_0(i,j)*dnh*Mh2o)))/dt;
        
        end
        
    
    elseif INDC2(i,j)==6 || INDC2(i,j)==7
        
       Rw=-Rw+vb(i,j)*((phi(i,j)*sw(i,j)*dnw(i,j)*(1-cl(i,j))+phi(i,j)*si(i,j)*dni+phi(i,j)*sh(i,j)*dnh*Mh2o)-...
            ((phi_0(i,j)*sw_0(i,j)*dnw_0(i,j)*(1-cl_0(i,j))+phi_0(i,j)*si_0(i,j)*dni+phi_0(i,j)*sh_0(i,j)*dnh*Mh2o)))/dt;
    end
end



if i==COL
    Rw=Rw-dy(i,j)*dz(i,j)*qw/(365*24*3600);
end


