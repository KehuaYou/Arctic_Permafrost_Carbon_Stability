
global COL ROW
global INDC1 INDC2 INDC3 
global pw
global  sh  sg  si
global sh_old si_old sg_old
global cl  cm
global  salinity
global  T
global p_salinity T_salinity
global go_back


for i=1:COL
    for j=1:ROW
        %*************************** smoothly change indexes***************
        if INDC2(i,j) == 2
            INDC2(i,j)=3;
            go_back=go_back+0;
        end
        
        if INDC2(i,j) == 4
            INDC2(i,j)=5;
            go_back=go_back+0;
        end
        
        
        if INDC2(i,j) == 6
            INDC2(i,j)=7;
            go_back=go_back+0;
        end
        
        if  INDC3(i,j)==3
            INDC3(i,j)=4;
        end
        
        
        %******************one phase changes to two phases****************
        %********  one phase can only be L here          *****************
        %********  two phases: L+G, L+H, L+I             *****************
        %*****************************************************************
        if INDC2(i,j)==1 && INDC3(i,j) == 1
            if INDC1(i,j) == 0
                if cm(i,j) > methane_solubility(pw(i,j)/1e6,T(i,j),2*cl(i,j)/0.11099/(1-cl(i,j)))*0.016    % L changes to L+G
                    INDC2(i,j)=2;
                    go_back=go_back+1;
                end
            else
                if cm(i,j) > hydrate_solubility(pw(i,j)/1e6,T(i,j),2*cl(i,j)/0.11099/(1-cl(i,j)))*0.016   % L changes to L+H
                    INDC2(i,j)=2;
                    go_back=go_back+1;
                end
            end
            
            if T(i,j)<=-cl(i,j)*(164.49*cl(i,j)+49.462)     % L changes to L+I
                INDC2(i,j)=2;
                INDC3(i,j)=3;
                go_back=go_back+1;
            end
        end
        %****************** two phases change to three phases*****************
        %******** two phases:  L+G, L+H, L+I            **********************
        %******** three phases: L+G+I, L+G+H, L+I+H, I+H+G *******************
        %*********************************************************************
        
        if INDC2(i,j)==3
            if INDC1(i,j)==0
                if INDC3(i,j)==1 && cl(i,j) < interpolation2(p_salinity,T_salinity,salinity,pw(i,j)/1e6,T(i,j))  % L+G changes to L+H+G
                    INDC2(i,j)=4;
                    go_back=go_back+1;
                end
                if INDC3(i,j)==2 && pw(i,j)/1e6 > interpolation_pressure(p_salinity,T_salinity,salinity,cl(i,j),T(i,j))  % I+G changes to I+H+G
                    INDC2(i,j)=4;
                    go_back=go_back+1;
                end
                if INDC3(i,j)==4 && cm(i,j) > methane_solubility(pw(i,j)/1e6,T(i,j),2*cl(i,j)/0.11099/(1-cl(i,j)))*0.016  % L+I changes to L+I+G
                    INDC2(i,j)=4;
                    go_back=go_back+1;
                end
            else
                if INDC3(i,j)==1 && cl(i,j) > interpolation2(p_salinity,T_salinity,salinity,pw(i,j)/1e6,T(i,j)) % L+H changes to L+H+G
                    INDC2(i,j)=4;
                    go_back=go_back+1;
                end
                if INDC3(i,j)==2 && pw(i,j)/1e6 < interpolation_pressure(p_salinity,T_salinity,salinity,cl(i,j),T(i,j))  % I+H changes to I+H+G
                    INDC2(i,j)=4;
                    go_back=go_back+1;
                end
                if  INDC3(i,j)==4 && cm(i,j) > hydrate_solubility(pw(i,j)/1e6,T(i,j),2*cl(i,j)/0.11099/(1-cl(i,j)))*0.016  % L+I changes to L+I+H
                    INDC2(i,j)=4;
                    go_back=go_back+1;
                end
            end
            
            
            if INDC3(i,j)==1 && T(i,j)<=-cl(i,j)*(164.49*cl(i,j)+49.462)  % 1) L+G changes to L+I+G; 2)L+H changes to L+I+H
                INDC2(i,j)=4;
                INDC3(i,j)=3;
                go_back=go_back+1;
            end
            
            if INDC3(i,j)==2 && T(i,j)>=-cl(i,j)*(164.49*cl(i,j)+49.462)  % 1) I+G changes to L+I+G; 2) I+H changes to L+I+H
                INDC2(i,j)=4;
                INDC3(i,j)=3;
                go_back=go_back+1;
            end
        end
        
        %****************** three phases change to four phases****************
        %******** three phases: L+G+I, L+G+H, L+I+H, I+H+G *******************
        %******** four phases: L+I+H+G                     *******************
        %*********************************************************************
        
        if INDC2(i,j)==5
            if INDC3(i,j)==1 && T(i,j)<=-cl(i,j)*(164.49*cl(i,j)+49.462) % L+H+G changes to L+I+H+G
                INDC2(i,j)=6;
                INDC3(i,j)=3;
                go_back=go_back+1;
            end
            
            if INDC3(i,j)==2 && T(i,j)>=-cl(i,j)*(164.49*cl(i,j)+49.462)    % I+H+G changes to L+I+H+G
                INDC2(i,j)=6;
                INDC3(i,j)=3;
                go_back=go_back+1;
            end
            
            if INDC3(i,j)==4
                if INDC1(i,j)==0 && cl(i,j) < interpolation2(p_salinity,T_salinity,salinity,pw(i,j)/1e6,T(i,j)) % L+I+G changes to L+I+H+G
                    INDC2(i,j)=6;
                    go_back=go_back+1;
                end
                
                if INDC1(i,j)==1 && cl(i,j) > interpolation2(p_salinity,T_salinity,salinity,pw(i,j)/1e6,T(i,j)) % L+I+H changes to L+I+H+G
                    INDC2(i,j)=6;
                    go_back=go_back+1;
                end
            end
        end
        
        %**************Phase number decreases*********************
        
        if sh(i,j)<0 && sh_old(i,j)>0       % hydrate disappear
            if INDC2(i,j)>5             % I+L+G+H  changes to I+L+G
                INDC2(i,j)=4;
                INDC1(i,j)=0;
            elseif INDC2(i,j)>3 % I+H+G  changes to I+G, L+H+G changes to L+G, L+I+H -> L+I
                INDC2(i,j)=2;
                INDC1(i,j)=0;
            else               % L+H -> H, I+H->I
                INDC2(i,j)=1;
            end
            go_back=go_back+1;
        end
        
        if sg(i,j)<0  && sg_old(i,j)>0     % gas disappear
            if INDC2(i,j)>5             % I+L+G+H changes to I+L+H
                INDC2(i,j)=4;
                INDC1(i,j)=1;
            elseif INDC2(i,j)>3        % I+H+G changes to I+H, L+H+G changes to L+H, L+I+G -> L+I
                INDC2(i,j)=2;
                INDC1(i,j)=1;
            else
                INDC2(i,j)=1;
            end
            go_back=go_back+1;
        end
        
        if si(i,j)<0 && si_old(i,j)>0   % ice disappear
            if INDC2(i,j)>5         % I+L+G+H changes to L+H+G
                INDC2(i,j)=4;
            elseif INDC2(i,j)>3     % I+L+H changes to L+H or I+L+G changes to L+G, I+H+G -> H+G
                INDC2(i,j)=2;
            else                        % I+L changes to L, I+G -> G, I+H->H
                INDC2(i,j)=1;
            end
            INDC3(i,j)=1;
            go_back=go_back+1;
        end
    end
end