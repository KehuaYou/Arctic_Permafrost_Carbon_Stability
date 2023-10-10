function y=interpolation_pressure(p_salinity,T_salinity,salinity,cl,T)
[row,col]=size(p_salinity);

if cl<=salinity(1)
    xarr=p_salinity(:,1);
elseif cl>=salinity(col)
    xarr=p_salinity(:,col);
else
    for i=1:col-1
        if cl>salinity(i) && cl<=salinity(i+1)
            xarr=(cl-salinity(i))/(salinity(i+1)-salinity(i))*(p_salinity(:,i+1)-p_salinity(:,i))+p_salinity(:,i); % linear interpolation between P and salinity
        end
    end
end


if T<=T_salinity(1)
    y=xarr(1);
elseif T>=T_salinity(row)
    y=xarr(row);
else
    for i=1:row-1
        if T>T_salinity(i) && T<=T_salinity(i+1)
            y=(T-T_salinity(i))/(T_salinity(i+1)-T_salinity(i))*(xarr(i+1)-xarr(i))+xarr(i);    % linear interpolation between P and T
        end
    end
end
 end
 
    
    
 

