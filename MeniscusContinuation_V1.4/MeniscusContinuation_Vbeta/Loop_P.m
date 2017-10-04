function [R,Z,P,V] = Loop_P(R,Z,P,V,dP,Npas)

global Bo RR dpdz alphac Weta Wxi Vref; 


conv = 1;
Pstart = P;
PP = Pstart;
Vtab = [];
Ptab = [];
itloop = 0;

while ((itloop<=Npas)&&(conv==1))
[R,Z,P,V,conv] = Newton_P(R,Z,P,V,PP);
if(conv)
Vtab = [Vtab V];
Ptab = [Ptab PP];
end

PP = PP+dP;
itloop = itloop+1;
end



figure(21);
hold on;
plot(Ptab,Vtab/Vref,'r');


figure(10);
hold on;
plot(R,Z,'r');

