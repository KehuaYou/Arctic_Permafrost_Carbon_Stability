function y=interp(xarr,yarr,NUM,x,flag)
if x <= xarr(1)
    y=yarr(1);
elseif x >= xarr(NUM)
    y=yarr(NUM);
else
    for i=1:NUM-1
        if x > xarr(i) && x <= xarr(i+1)
            if flag==0
                y=(log10(x)-log10(xarr(i)));
                y=y/(log10(xarr(i+1))-log10(yarr(i)));
                y=y*(log10(yarr(i+1))-log10(yarr(i)))+log10(yarr(i));
                y=10^y;
            elseif flag==1
                y=(x-xarr(i))/(xarr(i+1)-xarr(i))*(yarr(i+1)-yarr(i))+yarr(i);
            else
                'You need to specify flag'
            end
        end
    end
end