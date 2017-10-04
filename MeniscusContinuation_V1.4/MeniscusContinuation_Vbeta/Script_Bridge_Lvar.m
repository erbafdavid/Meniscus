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


% geometrie initiale : pont cylindrique
t = [0:.02:1];
R0 = RR*ones(size(t));

for L = 1:.2:8

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

 [alpha,ds,s0,Ka,Kb] = calcgeom(R,Z);
[eta,P1,RES,lambdaminV,lambdamaxV,lambdaminP,lambdamaxP] = calceta_V(R,Z,P,alpha,ds,s0,Ka,Kb,dpdz,V)
figure(32);
hold on;
plot(L,lambdaminP,'xr',L,real(lambdamaxP),'or',L,real(lambdaminV),'xb',L,real(lambdamaxV),'ob');

end




