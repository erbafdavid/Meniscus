


clear all;
close all;

% construction de la courbe P/V pour une goutte pendante,
% Bo = 1
% (adimensionalisation par la longueur capilaire)


global Bo dpdz alphac Weta Wxi verbosity beta discretization nbdim Vref;
global tabP tabV tab1; 


verbosity = 0;
discretization = 'FD';
beta = -1;
nbdim = 3;
Vref = 1;


dpdz = -1;
Weta = 1;
Wxi = 1;
t = [0:.01:1];


tabV = [0;0;0];
tabP = [0;0;0];
tab1 = [0;0;0];


for Bo = 1:.1:5.01

Bo
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = 0.2*Bo^2;

[R,Z,P,V] = Newton_P(R,Z,P,V,.01);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.005,100);
end

figure(71);
plot(tabP(1,2:length(tabP)),tabP(2,2:length(tabP)),'r');

figure(72);
plot(tabV(1,2:length(tabP)),tabV(3,2:length(tabP))./tabV(1,2:length(tabP)).^2/pi,'k--');
hold on;
plot(tab1(1,2:length(tabP)),tab1(3,2:length(tabP))./tab1(1,2:length(tabP)).^2/pi,'r');
plot(tabP(1,2:length(tabP)),tabP(3,2:length(tabP))./tabP(1,2:length(tabP)).^2/pi,'b:');



tabV = [0;0;0];
tabP = [0;0;0];
tab1 = [0;0;0];

for Bo = 1:-.05:.1

Bo
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = 0.1*Bo;

[R,Z,P,V] = Newton_P(R,Z,P,V,.0001);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV*Bo*2,12); % pour faible Bo
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.005,100);
end

figure(71);
hold on;

plot(tabP(1,2:size(tabP,2)),tabP(2,2:size(tabP,2)),'r');

%tabP = [tabP,[0;100;0]];


figure(72);
plot(tabV(1,2:length(tabV)),tabV(3,2:length(tabV))./tabV(1,2:length(tabV)).^2/pi,'k--');
hold on;
plot(tab1(1,2:length(tab1)),tab1(3,2:length(tab1))./tab1(1,2:length(tab1)).^2/pi,'r');
plot(tabP(1,2:size(tabP,2)),tabP(3,2:size(tabP,2))./tabP(1,2:size(tabP,2)).^2/pi,'b:');


figure(71);
xlabel('Bo');
ylabel('P');
box on;

figure(72);
xlabel('Bo');
ylabel('H');
box on;




