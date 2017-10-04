% construction de la courbe P/V pour un pont
% (adimensionalisation par le rayon)

% CAS L=2, SANS GRAVITE

global Bo RR dpdz alphac Weta Wxi verbosity;
Bo = 1e12;
RR = 1;
dpdz = 0;
Weta = 1;
Wxi = 0.2;
verbosity = 0;
close all;

L = 3;



% solution limite theorique : deux portions de sphere de rayon Rc tel que
% Rc^2 = R^2 + (Rc-L/2)^2 
F = @(x)(x^2-1-(x-L/2)^2);
Rc = fzero(F,1.5)
%  Rc = 1
Kl = 2/Rc;
Vl = 2*pi/3*(L/2)^2*(3*Rc-L/2);% formule volume d'une calotte spherique
figure(21);
plot(Kl,Vl,'go','MarkerSize',10);

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
plot(R,Z,'r');
figure(21);
plot(P,V,'ro','MarkerSize',10);

% continuation en diminuant P 
[R,Z,P,V] = Loop_P(R,Z,P,V,-.02,5);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.02,5);


% au passage on attrape la solution d'equilibre final 
% (volume = celui des deux portions de spheres)
[Rf,Zf,Pf,Vf] = Newton_V(R,Z,P,V,Vl);
figure(21);
plot(Pf,Vf,'go','MarkerSize',10);
figure(10);
plot(Rf,Zf,'g')

% loop V dans le virage
[R,Z,P,V] = Loop_V(R,Z,P,V,-.05,5);
[R,Z,P,V] = Loop_V(R,Z,P,V,-.05,20);

%remonte en P
[R,Z,P,V] = Loop_P(R,Z,P,V,+.03,30);
[R,Z,P,V] = Loop_P(R,Z,P,V,+.03,30);
[R,Z,P,V] = Loop_V(R,Z,P,V,+.1,50);


if(1==0)
[R,Z,P,V] = Loop_P(R,Z,P,V,+.1,10);
[R,Z,P,V] = Loop_P(R,Z,P,V,+.02,10);
[R,Z,P,V] = Loop_P(R,Z,P,V,+.004,10);

%solution limite deux spheres qui se touchent
[Rlim,Zlim,Plim,Vlim] = Newton_P(R,Z,P,V,P);
figure(10);
plot(Rlim,Zlim,'g')

% continuation en montant P 
[R,Z,P,V] = Loop_V(R0,Z0,1,0,+.3,5);
[R,Z,P,V] = Loop_V(R,Z,P,V,+.6,30);

end

% mise en forme des figures

figure(21);
grid on;
xlabel('P')
ylabel('V')

figure(10);
grid on;
xlabel('R')
ylabel('Z')



