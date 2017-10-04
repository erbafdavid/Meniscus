

clear all;
close all;

% test for pending drop with prescribed contact angle
global Bo RR dpdz alphac Weta Wxi verbosity beta discretization nbdim;


nbdim = 3; % 3d or 2d
Bo = 1;
RR = 1;
dpdz = .0;
Weta = 1;
Wxi = 1;

beta = pi/180*72;
betainit = beta;
discretization = 'FD';
verbosity = 0;


RRR = .2;
t = [0:.02:1];
R0 = RRR*sin(betainit*t);
Z0 = -RRR*(cos(betainit*t)-cos(betainit));


[alpha,ds,s0,Ka,Kb] = calcgeom(R0,Z0);
alpha(length(alpha))

P = 2./RRR;
V = 0;



dpdz = 1;
[R,Z,P,V,conv] = Newton_P(R0,Z0,P,V,P);
P
V

[R,Z,P,V] = Loop_P(R,Z,P,V,-.4,10);

[R,Z,P,V] = Loop_V(R,Z,P,V,.2,100);

[R,Z,P,V] = Loop_V(R,Z,P,V,.01,20);

[R,Z,P,V] = Loop_P(R,Z,P,V,-.002,20);

[R,Z,P,V] = Loop_V(R,Z,P,V,-.005,30);

[R,Z,P,V] = Loop_V(R,Z,P,V,-.05,100);









