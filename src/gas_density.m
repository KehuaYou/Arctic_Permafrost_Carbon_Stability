function density=gas_density(P,T)
R=83.1441;
P=P*10;
T=273.15+T;
% P is in bar and T is in K. 
Tc=190.6;
Pc=46.41;
Pr=P/Pc;
Tr=T/Tc;
Vr=calc_Vr(Pr,Tr);
Z=Pr*Vr/Tr;
tmp_v=Z*R*T/(P*1e5)/10;
density=0.016/tmp_v;



    



