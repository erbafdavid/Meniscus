


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

Bo = 1;
RR = Bo;
dpdz = 1;
Weta = 1;
Wxi = 1;



t = [0:.05:1];
R0 = Bo*t;
Z0 = 0*t;



R = R0;
Z = Z0;
P = .01;
V = 0;


%[R,Z,P,V] = Newton_P(R0,Z0,P,V,P);

%if(1==0)
[R,Z,P,V] = Loop_P(R,Z,P,V,.1,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,.05,10);
[R,Z,P,V] = Loop_V(R,Z,P,V,.1,10);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.05,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-.05,100);




if (1==0)
    
R = R0;
Z = Z0;
P = .01;
V = 0;

discretization = 'FE';

%[R,Z,P,V] = Newton_P(R0,Z0,P,V,P);
%if(1==0)
[R,Z,P,V] = Loop_P(R,Z,P,V,.1,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,.05,10);
[R,Z,P,V] = Loop_V(R,Z,P,V,.1,10);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.05,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-.05,100);

end

