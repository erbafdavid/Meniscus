


clear all;
close all;

% construction de la courbe P/V pour une goutte pendante,
% Bo = 1
% (adimensionalisation par la longueur capilaire)


global RR dpdz alphac Weta Wxi verbosity beta discretization nbdim Vref;


verbosity = 0;

beta = -1;
nbdim = 3;
Vref = 1;

Bo = 1;
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
%%%
Bo = 1;
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = .1;


discretization = 'FD';

[R,Z,P,V] = Loop_V(R,Z,P,V,1.399/20,20);


[R,Z,P,V] = Newton_P(R,Z,P,V,.001);
[R,Z,P,V] = Loop_P(R,Z,P,V,.05,10);

[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.05,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-dV,20);
end




    
R = R0;
Z = Z0;
P = .01;
V = 0;

discretization = 'FD';


%[R,Z,P,V] = Newton_P(R,Z,P,V,.0);
[R,Z,P,V] = Loop_V(R,Z,P,V,1.399/20,20);

%[R,Z,P,V] = Loop_P(R,Z,P,V,-.05,100);

%[R,Z,P,V] = Loop_V(R,Z,P,V,-dV,20);
%[R,Z,P,V] = Loop_P(R,Z,P,V,.05,100);
%[R,Z,P,V] = Newton_P(R,Z,P,V,V);
