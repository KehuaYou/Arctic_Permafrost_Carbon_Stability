function fugacity=calc_fugacity(P0,T0,aw)
x0=0.5*P0*10;
x=x0;
loop=1;
iteration=0;
eps=1e-4;
while loop==1
    iteration=iteration+1;
    y0=fugacity_methane(x,P0,T0,aw);
    B=-y0;
    x=x+eps;
    y1=fugacity_methane(x,P0,T0,aw);
    A=(y1-y0)/eps;
    delta=A\B;
    x=x+delta;
    if iteration > 20
        loop=0;
    end
end
fugacity=x;
    
    
   
    
    
    
