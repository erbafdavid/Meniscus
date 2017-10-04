% construction de la courbe P/V pour un pont
% (adimensionalisation par le rayon)

% CAS L=8, SANS GRAVITE

close all;

global Bo RR dpdz alphac Weta Wxi verbosity beta discretization nbdim Vref;


nbdim = 3; % 3d or 2d
Bo = 1;
RR = 1;
dpdz = .0;
Weta = 1;
Wxi = 1;

beta = -1;
discretization = 'FD';
verbosity = 0;

L = 8;
Vref = pi*L;




% solution limite theorique : deux portions de sphere de rayon Rc tel que
% Rc^2 = R^2 + (Rc-L/2)^2 
F = @(x)(x^2-1-(x-L/2)^2);
Rc = fzero(F,1.5)
%  Rc = 
Kl = 2/Rc;
Vl = 2*pi/3*(L/2)^2*(3*Rc-L/2);% formule volume d'une calotte spherique
figure(21);
plot(Kl,Vl/Vref,'go','MarkerSize',10);




% geometrie initiale : pont cylindrique
t = [0:.005:1];
N = length(t);
R0 = RR*ones(size(t));
Z0 = L*t;

R = R0;
Z = Z0;
P = 1;
V = 1; 
[R,Z,P,V] = Newton_P(R0,Z0,P,V,1);
figure(10);
hold on;
plot(R,Z,'m');
plot(0*R,Z,'k-.');
plot(1,0,'k.');
plot(1,L,'k.');


figure(21);
plot(P,V/Vref,'ro','MarkerSize',10);




% continuation en diminuant P 
[R,Z,P,V] = Loop_P(R,Z,P,V,+.01,17);

%[Raut,Zaut,Paut,Vaut] = Newton_P(R,Ztrans,P,V,P); % marche plus ???

figure(31);
grid on;
% Au passage on choppe un point sur la branche symetrique à deux cols

[Rxx,Zxx,Pxx,Vxx] = Newton_P(R,Z,P,V,P); 
[Raut,Zaut,Paut,Vaut] = Newton_P(R,Z,P,V,P); 
Raut = [Raut(round(N/3)+1:N) , Raut(round(N/3)+1:2*round(N/3))];
[Raut,Zaut,Paut,Vaut] = Newton_P(Raut,Zaut,Paut,Vaut,Paut);
Raut = Raut+1-Raut(1);
[Raut,Zaut,Paut,Vaut] = Newton_P(Raut,Zaut,Paut,Vaut,Paut);

if(1==1)

% virage
%[R,Z,P,V] = Loop_P(R,Z,P,V,+.02,30);
%remontee
[R,Z,P,V] = Loop_V(Rxx,Zxx,Pxx,Vxx,1,10);
[R,Z,P,V] = Loop_V(R,Z,P,V,2,10);
[R,Z,P,V] = Loop_V(R,Z,P,V,1,30);



%pause;


%solution limite deux spheres qui se touchent
[Rlim,Zlim,Plim,Vlim] = Newton_P(R,Z,P,V,P);
figure(10);
plot(Rlim,Zlim,'g')





% au passage on attrape le point de bifurcation
Pbif = 0.7775;
[Rbif,Zbif,Pbif,Vbif] = Newton_P(R0,Z0,1,1,Pbif);
figure(21);
plot(Pbif,Vbif/Vref,'bo','MarkerSize',10);
figure(10);
plot(Rbif,Zbif,'b')

% au passage on attrape la solution d'equilibre final 
% (volume = celui des deux portions de spheres)
[Rf,Zf,Pf,Vf] = Newton_V(Rbif,Zbif,Pbif,Vbif,Vl);
figure(21);
plot(Pf,Vf/Vref,'go','MarkerSize',10);
figure(10);
plot(Rf,Zf,'g')


% continuation en montant V  pour les ponts bombes
[R,Z,P,V] = Loop_P(R0,Z0,1,0,-.02,10);
[R,Z,P,V] = Loop_V(R,Z,P,V,(7*Vref-V)*.05,20);




% branche asymetrique : point limite
alpha = asin(2*RR/L);
tt = pi*t-alpha;
Ra = L/2*abs(sin(tt));
Za = L/2*(((tt>0).*cos(tt)+(tt<0).*(2-cos(tt)))+cos(alpha));
Va = L^3*pi/6;
Pa = 4/L;

figure(21);
hold on;
plot(Pa,Va/Vref,'bo','MarkerSize',10);

figure(10);
plot(Ra,Za,'b');
hold on;


% pour finir : branche asymetrique
dpdz =-0.001;
[Rg,Zg,Pg,Vg] = Loop_P(Rbif,Zbif,Pbif,Vbif,-.02,10);
%[Rx,Zx,Px,Vx] = Newton_P(R,Z,P,V,P); % pour attrapper la branche asymetrique au retour



dpdz =0;
[R,Z,P,V] = Loop_P(Rg,Zg,Pg,Vg,(Pbif-Pg)/20,19);



[R,Z,P,V] = Loop_V(Rg,Zg,Pg,Vg,(Va-Vg)*.05,19);

end

% Branche à deux cols 

figure(11);
plot(Raut,Zaut);

Pbif2 = 1.48;

[R,Z,P,V] = Loop_P(Raut,Zaut,Paut,Vaut,(Pbif2-P)/10,10);
[Rbif2,Zbif2,Pbif2,Vbif2] = Newton_P(R,Z,P,V,P);
figure(21);
hold on;
plot(Pbif2,Vbif2/Vref,'co','MarkerSize',10);
%vers la droite
[R,Z,P,V] = Loop_P(R,Z,P,V,.02,30);
[R,Z,P,V] = Loop_V(R,Z,P,V,+.05*Vref,30);
% vers le haut
[R,Z,P,V] = Loop_P(Raut,Zaut,Paut,Vaut,-.02,30);
[R,Z,P,V] = Loop_V(R,Z,P,V,+.2*Vref,30);

% branche asymetrique
dpdz =-0.001;
[Rg2,Zg2,Pg2,Vg2] = Loop_P(Rbif2,Zbif2,Pbif2,Vbif2,-.02,10);
dpdz = 0;
[R,Z,P,V] = Loop_V(Rg2,Zg2,Pg2,Vg2,.05*Vref,50);


%points de termination
F = @(x)(4*x-2*sqrt(x^2-1)-L);
Rc = fzero(F,3);
Kl2a = 2/Rc;
Vl2a = 4*pi/3*Rc^3+2*pi/3*(L/2-Rc)^2*(3*Rc-(L/2-Rc));% formule volume d'une calotte spherique
figure(21);
plot(Kl2a,Vl2a/Vref,'mo','MarkerSize',10);

%points de termination
F = @(x)(4*x+2*sqrt(x^2-1)-L);
Rc = fzero(F,3);
Kl2b = 2/Rc;
Vl2b = 4*pi/3*Rc^3+2*pi/3*(L/2-Rc)^2*(3*Rc-(L/2-Rc));% formule volume d'une calotte spherique
figure(21);
plot(Kl2b,Vl2b/Vref,'mo','MarkerSize',10);

Rc = L/4;
Kl2c = 2/Rc;
Vl2c = 4*pi/3*Rc^3*2;
figure(21);
plot(Kl2c,Vl2c/Vref,'co','MarkerSize',10);

% mise en forme des figures

figure(21);
grid on;
xlabel('P')
ylabel('V')

figure(10);
grid on;
xlabel('R')
ylabel('Z')



