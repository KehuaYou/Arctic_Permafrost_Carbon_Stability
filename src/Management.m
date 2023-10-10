global ROW COL
global sg sh si sw cm
global INDC2
global t_flag cycle

for i=1:COL
    for j=1:ROW
        if sg(i,1)<-0.0001
            sg(i,1)=0;
        end
        
        if sh(i,1)<-0.0001
            sh(i,1)=0;
        end
        
        if si(i,1)<-0.0001
            si(i,1)=0;
        end
        
        if sw(i,1)>1
            sw(i,1)=1;
        end
        
        if sg(i,1)<1e-6 && sg(i,1)>0 && INDC2(i,1)==3
            sg(i,1)=0;
            sw(i,1)=1;
            INDC2(i,1)=1;
        end
        if cm(i,1)<0
            cm(i,1)=0;
        end
        
        if sw(i,1)<0
            t_flag=1;
        end
        
        if cycle>4
            t_flag=1;
        end
    end
end