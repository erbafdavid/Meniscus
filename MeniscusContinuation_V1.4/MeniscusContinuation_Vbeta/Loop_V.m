function [R,Z,P,V] = Loop_V(R,Z,P,V,dV,Npas)

global Bo RR dpdz alphac Weta Wxi verbosity Vref; 

conv = 1;
VV = V;
Vtab = [];
Ptab = [];
itloop = 0;

while ((itloop<=Npas)&&(conv==1))

[R,Z,P,V,conv] = Newton_V(R,Z,P,V,VV);
if(conv)
Vtab = [Vtab VV];
Ptab = [Ptab P];
end;

VV = VV+dV;
itloop = itloop+1;
end

figure(21);
hold on;
plot(Ptab,Vtab/Vref,'r');



figure(10);
hold on;
plot(R,Z,'r');