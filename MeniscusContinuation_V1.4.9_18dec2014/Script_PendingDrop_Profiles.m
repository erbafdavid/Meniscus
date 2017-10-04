


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
dpdz = 1;
Weta = 1;
Wxi = 1;
t = [0:.01:1];

GG = 1;

gcf = figure(10);

set(gcf,'PaperUnits','centimeters')
xSize = 5.5; ySize = 5.9;
xLeft = 0; yTop = 0;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[1 1 xSize*100 ySize*100])

width = 1/xSize; height = 3/ySize;
x = .1/xSize; y = 1-(3+.1)/ySize; 
axes('position',[x y width height])
plot([0 0],[-2.8,.2],'-.k')
axis([0 1. -2.8 .2]);
ax = gca;
ax.XColor = 'w';
ax.YColor = 'w';

figure(21);
subplot(1,2,1);
box on;
xlabel('P');
ylabel('V');

%%%
%Bo = .2;
%R0 = Bo*t;
%Z0 = 0*t;
%R = R0;
%Z = Z0;
%P = .01;
%V = 0;
%dV = 0.02;

%[R,Z,P,V] = Newton_P(R,Z,P,V,.01);
%[R,Z,P,V] = Loop_V(R,Z,P,V,dV/10,50);
%[R,Z,P,V] = Loop_V(R,Z,P,V,2*dV,100);
%[R,Z,P,V] = Loop_P(R,Z,P,V,-.005,100);
%[R,Z,P,V] = Loop_V(R,Z,P,V,-dV/5,10);



%%%
Bo = .4;
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = 0.1;

[R,Z,P,V] = Newton_P(R,Z,P,V,.001);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV/5,10);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.01,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-dV,20);

%%%
Bo = .6;
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = 0.06;

[R,Z,P,V] = Newton_P(R,Z,P,V,.001);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.01,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-dV,20);


%%%
Bo = 1;
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = .1;

[R,Z,P,V] = Newton_P(R,Z,P,V,.001);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.01,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-dV,20);



figure(10);
width = 3.8/xSize; height = 3/ySize;
x = (.1+1+.5)/xSize; y = 1-(3+.1)/ySize; 
axes('position',[x y width height])
plot([0 0],[-2.8,.2],'-.k')
axis([0 3.8 -2.8 .2])
ax = gca;
ax.XColor = 'w';
ax.YColor = 'w';


%%%
Bo = 1.6;
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = .2;

[R,Z,P,V] = Newton_P(R,Z,P,V,.001);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.01,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-dV,30);

figure(21);
subplot(1,2,2);
box on;
xlabel('P');
ylabel('V');


%%%
Bo = 2.4;
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = .2;

[R,Z,P,V] = Newton_P(R,Z,P,V,.001);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.01,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-dV,30);


%%%
Bo = 3.5;
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = .2;

[R,Z,P,V] = Newton_P(R,Z,P,V,.001);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.01,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-dV,30);

figure(10);
width = 5/xSize; height = 2.2/ySize;
x = .1/xSize; y = 0.1/ySize; 
axes('position',[x y width height])
plot([0 0],[-2,.2],'-.k')
axis([0 5 -2 .2])
ax = gca;
ax.XColor = 'w';
ax.YColor = 'w';

%%%
Bo = 3.8;
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = .2;

[R,Z,P,V] = Newton_P(R,Z,P,V,.001);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.01,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-dV,30);



%%%
Bo = 4.5;
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = .05;

[R,Z,P,V] = Newton_V(R,Z,P,V,-.001);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.002,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-dV,30);

%%%
Bo = 5;
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = .05;

[R,Z,P,V] = Newton_V(R,Z,P,V,-.001);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.002,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-dV,30);




