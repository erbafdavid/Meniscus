


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
dpdz = -1;
Weta = 1;
Wxi = 1;
t = [0:.01:1];

GG = 1;

gcf = figure(10);

set(gcf,'PaperUnits','centimeters')
xSize = 14; ySize = 7.5;
xLeft = 0; yTop = 0;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[1 1 xSize*60 ySize*60])

width = 1.5/xSize; height = 3/ySize;
x = .35/xSize; y = 1-(.5+3)/ySize; 
axes('position',[x y width height])
plot([0 0],[-.5,2.5],'-.k')
axis([0 1.5 -.5 2.5])
ax = gca;
ax.XColor = 'w';
ax.YColor = 'w';

figure(21);
subplot(2,2,1);
box on;
    
%%%
Bo = .2;
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = 0.02;

[R,Z,P,V] = Newton_P(R,Z,P,V,.01);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV/10,50);
[R,Z,P,V] = Loop_V(R,Z,P,V,2*dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.005,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-dV/5,10);

%%%
Bo = .4;
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = 0.1;

[R,Z,P,V] = Newton_P(R,Z,P,V,.01);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV/5,10);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.05,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-dV,10);

%%%
Bo = .6;
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = 0.1;

[R,Z,P,V] = Newton_P(R,Z,P,V,.01);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.05,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-dV,10);

figure(10);
width = 4/xSize; height = 3/ySize;
x = (.35+1.5+.5)/xSize; y =1-(.5+3)/ySize; 
axes('position',[x y width height])
plot([0 0],[-.5,2.5],'-.k')
axis([0 4 -.5 2.5])
ax = gca;
ax.XColor = 'w';
ax.YColor = 'w';

figure(21);
subplot(2,2,2);
box on;

%%%
Bo = 1;
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = Bo;

[R,Z,P,V] = Newton_P(R,Z,P,V,.01);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.05,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-dV,10);


%%%
Bo = 2;
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = Bo;

[R,Z,P,V] = Newton_P(R,Z,P,V,.01);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.05,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-dV,10);



%%%
Bo = 3;
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = Bo;

[R,Z,P,V] = Newton_P(R,Z,P,V,.01);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV/2,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.05,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-dV,10);





figure(10);
width = 7/xSize; height = 3/ySize;
x = (.35+1.5+.5+4+.5)/xSize; y =1-(.5+3)/ySize; 
axes('position',[x y width height])
plot([0 0],[-.5,2.5],'-.k')
axis([0 7 -.5 2.5])
ax = gca;
ax.XColor = 'w';
ax.YColor = 'w';

figure(21);
subplot(2,2,3);
box on;

%%%
Bo = 4;
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = Bo;

[R,Z,P,V] = Newton_P(R,Z,P,V,.01);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.05,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-dV,10);    
    
%%%
Bo = 5;
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = Bo;

[R,Z,P,V] = Newton_P(R,Z,P,V,.01);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.05,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-dV,10);

%%%
Bo = 6;
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = Bo;

[R,Z,P,V] = Newton_P(R,Z,P,V,.01);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.05,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-dV,10);

figure(10);
width = 13/xSize; height = 3/ySize;

x = .5/xSize; y =.5/ySize; 
axes('position',[x y width height])
plot([0 0],[-.5,2.5],'-.k')
axis([0 13 -.5 2.5])
ax = gca;
ax.XColor = 'w';
ax.YColor = 'w';

figure(21);
subplot(2,2,4);
box on;

%%%
Bo = 8;
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = Bo;

[R,Z,P,V] = Newton_P(R,Z,P,V,.01);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.05,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-dV,10);
%%%
Bo = 10;
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = Bo;

[R,Z,P,V] = Newton_P(R,Z,P,V,.01);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.05,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-dV,10);

%%%
Bo = 12;
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = Bo;

[R,Z,P,V] = Newton_P(R,Z,P,V,.01);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.05,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,-dV,10);


