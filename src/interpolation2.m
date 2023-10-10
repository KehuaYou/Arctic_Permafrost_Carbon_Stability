function y=interpolation2(p_salinity,T_salinity,salinity,p,T)
[row,col]=size(p_salinity);
if T <= T_salinity(1)
   xarr=p_salinity(1,:);
elseif T >= T_salinity(row)
   xarr=p_salinity(row,:);
else
    for i=1:row-1
        if T > T_salinity(i) && T <= T_salinity(i+1)
            xarr=(T-T_salinity(i))/(T_salinity(i+1)-T_salinity(i))*(p_salinity(i+1,:)-p_salinity(i,:))+p_salinity(i,:); %linear interpolation of pressure and temperature
        end
    end
end
 
 if p <= xarr(1)
    y=salinity(1);
 elseif p >= xarr(col)
    y=salinity(col);
 else
    for i=1:col-1
       if p > xarr(i) && p <= xarr(i+1)
          y=(p-xarr(i))/(xarr(i+1)-xarr(i))*(salinity(i+1)-salinity(i))+salinity(i); % linear interpolation of pressure and salinity
       end
    end
 end
 
    
    
 

