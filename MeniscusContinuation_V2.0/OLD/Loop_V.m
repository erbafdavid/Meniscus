function [R,Z,P,V] = Loop_V(R,Z,P,V,dV,Npas)

global Bo RR dpdz alphac Weta Wxi verbosity Vref tabV tabP tab1; 

conv = 1;
VV = V;
Vtab = [];
Ptab = [];
itloop = 0;
[R,Z,P,V,conv,lambdaminV,lambdaminP,lambdamin1] = Newton_V(R,Z,P,V,VV);


while ((itloop<=Npas)&&(conv==1))
lambdaminVans = lambdaminV;
lambdaminPans = lambdaminP;
lambdamin1ans = lambdamin1;
Vans = V;
Pans = P;

[R,Z,P,V,conv,lambdaminV,lambdaminP,lambdamin1] = Newton_V(R,Z,P,V,VV);
if(conv)
Vtab = [Vtab VV];
Ptab = [Ptab P];



if (lambdaminV*lambdaminVans<0&&abs(lambdaminV)~=1&&abs(lambdaminVans)~=1) % detection bifurcation V
    VbV = (V*lambdaminVans - Vans*lambdaminV)/(lambdaminVans-lambdaminV);
    PbV = (P*lambdaminVans - Pans*lambdaminV)/(lambdaminVans-lambdaminV);
    figure(21);
    hold on;
    plot(PbV,VbV/Vref,'+r','MarkerSize',10); 
    figure(10);
    hold on;
    plot(R,Z*sign(dpdz),'b')
    disp('Bifurcation V detectee !')
    PbV
    VbV
 %   if (VbV>.4) % bidouille faible Bo
    tabV = [tabV,[Bo ; PbV ; VbV]];
 %   else
%    tabP = [tabP,[Bo ; PbV ; VbV]];    
%    end %end bidouille
end

if (lambdaminP*lambdaminPans<0&&abs(lambdaminP)~=1&&abs(lambdaminPans)~=1) % detection bifurcation V
    VbP = (V*lambdaminPans - Vans*lambdaminP)/(lambdaminPans-lambdaminP);
    PbP = (P*lambdaminPans - Pans*lambdaminP)/(lambdaminPans-lambdaminP);
    figure(21);
    hold on;
    plot(PbP,VbP/Vref,'rx','MarkerSize',10); 
    figure(10);
    hold on;
    plot(R,Z*sign(dpdz),'k:')
    disp('Bifurcation P detectee !')
    PbP
    VbP
 
    [Bo ; PbP ; VbP]
 %   conv=0 pour arret lors du trace des courbes (L,V)
    
    
    tabP = [tabP,[Bo ; PbP ; VbP]];
end

if (lambdamin1*lambdamin1ans<0&&abs(lambdamin1)~=1&&abs(lambdamin1ans)~=1) % detection bifurcation V
    Vb1 = (V*lambdamin1ans - Vans*lambdamin1)/(lambdamin1ans-lambdamin1);
    Pb1 = (P*lambdamin1ans - Pans*lambdamin1)/(lambdamin1ans-lambdamin1);
    figure(21);
    hold on;
    plot(Pb1,Vb1/Vref,'*k','MarkerSize',10); 
    figure(10);
    hold on;
    plot(R,Z*sign(dpdz),'m')
    disp('Bifurcation m=1 detectee !')
    Pb1
    Vb1
    tab1 = [tab1,[Bo ; Pb1 ; Vb1]];
end


end;

VV = VV+dV;
itloop = itloop+1;
end

figure(21);
hold on;
plot(Ptab,Vtab/Vref,'r');



%figure(10);
%hold on;
%plot(R,Z,'r');