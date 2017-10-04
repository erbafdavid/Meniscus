% construction de la courbe P/V pour un pont
% (adimensionalisation par le rayon)

% CAS L=6, SANS GRAVITE

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

L = 6;
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
[R,Z,P,V] = Loop_P(R,Z,P,V,+.02,5);

figure(31);
grid on;


[R,Z,P,V] = Loop_P(R,Z,P,V,+.02,5);

[Rx,Zx,Px,Vx] = Newton_P(R,Z,P,V,P); % pour attrapper la branche asymetrique au retour

% virage
[R,Z,P,V] = Loop_P(R,Z,P,V,+.02,6);
%remontee
[R,Z,P,V] = Loop_V(R,Z,P,V,1,10);
[R,Z,P,V] = Loop_V(R,Z,P,V,1,10);
[R,Z,P,V] = Loop_V(R,Z,P,V,1,10);

%solution limite deux spheres qui se touchent
[Rlim,Zlim,Plim,Vlim] = Newton_P(R,Z,P,V,P);
figure(10);
plot(Rlim,Zlim,'g')



% au passage on attrape la solution d'equilibre final 
% (volume = celui des deux portions de spheres)
[Rf,Zf,Pf,Vf] = Newton_V(R0,Z0,1,1,Vl);
figure(21);
plot(Pf,Vf/Vref,'go','MarkerSize',10);
figure(10);
plot(Rf,Zf,'g')

% au passage on attrape le point de bifurcation
Pbif = 1.05;
[Rbif,Zbif,Pbif,Vbif] = Newton_P(R0,Z0,1,1,Pbif);
figure(21);
plot(Pbif,Vbif/Vref,'bo','MarkerSize',10);
figure(10);
plot(Rbif,Zbif,'b')


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
[Rg,Zg,Pg,Vg] = Loop_P(Rx,Zx,Px,Vx,-.02,15);

dpdz =0;
[R,Z,P,V] = Loop_P(Rg,Zg,Pg,Vg,(Pbif-Pg)/20,20);
figure(10);
plot(Rbif,Zbif,'b')

[R,Z,P,V] = Loop_V(Rg,Zg,Pg,Vg,(Va-Vg)*.05,19);




% mise en forme des figures

figure(21);
grid on;
xlabel('P')
ylabel('V')

figure(10);
grid on;
xlabel('R')
ylabel('Z')



