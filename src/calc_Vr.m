function Vr=calc_Vr(Pr,Tr)
x0=Tr/Pr;
x=x0;
loop=1;
iteration=0;
eps=1e-4;
while loop==1
    iteration=iteration+1;
    y0=duan(x,Pr,Tr);
    B=-y0;
    x=x+eps;
    y1=duan(x,Pr,Tr);
    A=(y1-y0)/eps;
    delta=A\B;
    x=x+delta;
    if iteration > 20
        loop=0;
    end
end
Vr=x;