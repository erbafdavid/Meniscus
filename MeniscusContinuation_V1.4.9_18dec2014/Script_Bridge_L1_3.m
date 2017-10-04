% construction de la courbe P/V pour un pont
% (adimensionalisation par le rayon)

% COMPARAISON AVEC LA THEORIE DANS LE CAS L=1.3, SANS GRAVITE
% cas dans lequel il existe deux solutions caténoidales à courbure nulle

close all;clear all;

global Bo RR dpdz alphac Weta Wxi beta discretization nbdim Vref verbosity;
nbdim = 3; % 3d or 2d
Bo = 1;
RR = 1;
dpdz = .0;
Weta = 1;
Wxi = 1;
beta = -1;
discretization = 'FD';
verbosity = 5;


L = 1.3;

Vref = pi*L;




% SOLUTIONS THEORIQUES CATENOIDALES :
% ces solutions sont de la forme r = a*cosh(z'/a) avec z'=z-L/2 et 
% ou a est la solution de a cosh(L/2a) = 1.

F = @(x)(x*cosh(.65/x)-1);
a1 = fzero(F,.64);
a2 = fzero(F,.46);
%a1 =   0.641607601777805;
%a2 =   0.461688844828455;

% solution limite theorique : deux portions de sphere de rayon Rc tel que
% Rc^2 = R^2 + (Rc-L/2)^2 
F = @(x)(x^2-1-(x-.65)^2);
Rc = fzero(F,1.5);
%  Rc = 1.094230769230769;
Kl = 2/Rc;
Vl = 2*pi/3*(L/2)^2*(3*Rc-L/2);% formule volume d'une calotte spherique


% geometrie initiale : pont cylindrique
t = [0:.005:1];
R0 = RR*ones(size(t));
Z0 = L*t;
figure(10);
plot(R0,Z0,'b')
hold on;
plot(R0*0,Z0,'k-.')
plot(1,0,'k.')
plot(1,L,'k.')



R = R0;
Z = Z0;
P = 1;
V = 1; 
[R,Z,P,V] = Newton_P(R0,Z0,P,V,1);
figure(10);
hold on;
plot(R,Z,'r');
figure(21);
plot(P,V/Vref,'rs','MarkerSize',10);

% continuation en diminuant P 
[R,Z,P,V] = Loop_P(R,Z,P,V,-.1,5);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.1,5);

% comparaison theorie premiere solution catenoidale
[Rc1,Zc1,Pc1,Vc1] = Newton_P(R,Z,P,V,0);
Rc1T = a1*cosh((Z0-L/2)/a1);
figure(11); 
plot(Rc1T,Z0,'b+');
figure(21);
plot(Pc1,Vc1/Vref,'bp','MarkerSize',10);
figure(10);
plot(Rc1,Zc1,'b')

% au passage on attrape la solution d'equilibre final 
% (volume = celui des deux portions de spheres)
[Rf,Zf,Pf,Vf] = Newton_V(R,Z,P,V,Vc1);
figure(21);
plot(Pf,Vf/Vref,'go','MarkerSize',10);
figure(10);
plot(Rf,Zf,'g')

% continuation en diminuant V 
[R,Z,P,V] = Loop_V(R,Z,P,V,-.05,6);
[R,Z,P,V] = Loop_V(R,Z,P,V,-.05,6);

% seconde solution catenoidale
[Rc2,Zc2,Pc2,Vc2] = Newton_P(R,Z,P,V,0);
% comparaison theorie
Rc2T = a2*cosh((Z0-L/2)/a2);
figure(10);
plot(Rc2,Zc2,'b')
figure(11); 
plot(Rc2T,Z0,'b+');
figure(21);
plot(Pc2,Vc2/Vref,'bp','MarkerSize',10);


% continuation en remontant P 
[R,Z,P,V] = Loop_P(R,Z,P,V,.1,10);
[R,Z,P,V] = Loop_P(R,Z,P,V,.1,10);
[R,Z,P,V] = Loop_P(R,Z,P,V,.02,20);
[Rlim,Zlim,Plim,Vlim] = Newton_P(R,Z,P,V,P);
figure(10);
plot(Rlim,Zlim,'g')

figure(21);
plot(Kl,Vl/Vref,'go','MarkerSize',10);

% mise en forme des figures

figure(21);
grid on;
xlabel('P')
ylabel('V')

figure(10);
grid on;
xlabel('R')
ylabel('Z')

% continuation en montant P 
[R,Z,P,V] = Loop_P(R0,Z0,1,0,+.1,20);
[R,Z,P,V] = Loop_V(R,Z,P,V,+.1,40);
[R,Z,P,V] = Loop_V(R,Z,P,V,+.02,40);
[R,Z,P,V] = Loop_V(R,Z,P,V,+.1,10);
%[R,Z,P,V] = Loop_V(R,Z,P,V,+.01,40);
