


clear all;
close all;

% construction de la courbe P/V pour une goutte pendante,
% Bo = 1
% (adimensionalisation par la longueur capilaire)

global Bo RR dpdz alphac Weta Wxi verbosity beta discretization nbdim Vref;


beta = -1; %contact angle (pinned if negative)



RR = 1;
dpdz = 0; % Bond = infinity  ; length scale = radius
Weta = 1;
Wxi = 1;
verbosity = 0;
nbdim = 3;
discretization = 'FD';
Vref = 1;


t = [0:.05:1];
R0 = RR*t;
Z0 = 0*t;

R = R0;
Z = Z0;
P = .01;
V = 0;


discretization = 'FD' ; % FE (finite elements) or FD (finite differences)
[R,Z,P,V] = Newton_P(R0,Z0,P,V,P);

%if(1==0)
[R,Z,P,V] = Loop_P(R,Z,P,V,.1,19);
[R,Z,P,V] = Loop_V(R,Z,P,V,(2*pi/3-V)/10.,30);


discretization = 'FE' ; % FE (finite elements) or FD (finite differences)
[R,Z,P,V] = Newton_P(R0,Z0,P,V,.01);
%R = R0;
%Z = Z0;
%P = .01;
%V = 0;
[R,Z,P,V] = Loop_P(R,Z,P,V,.1,19);
[R,Z,P,V] = Loop_V(R,Z,P,V,(2*pi/3-V)/10.,30);


