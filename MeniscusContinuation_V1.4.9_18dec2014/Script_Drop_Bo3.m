


clear all;
close all;

% construction de la courbe P/V pour une goutte pendante,
% Bo = 1
% (adimensionalisation par la longueur capilaire)


global RR dpdz alphac Weta Wxi verbosity beta discretization nbdim Vref;


verbosity = 0;
discretization = 'FD';
beta = -1;
nbdim = 3;
Vref = 1;

Bo = 3.6;
RR = Bo;
dpdz = 1;
Weta = 1;
Wxi = 1;



t = [0:.01:1];
R0 = Bo*t;
Z0 = 0*t;



R = R0;
Z = Z0;
P = .01;
V = 0;


%[R,Z,P,V] = Newton_P(R0,Z0,P,V,P);

%if(1==0)
[R,Z,P,V] = Newton_V(R,Z,P,V,.4);
[R,Z,P,V] = Loop_V(R,Z,P,V,.4,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.02,20);
[R,Z,P,V] = Loop_V(R,Z,P,V,-.2,20);

figure(31);
grid on;
figure (32);
grid on;

