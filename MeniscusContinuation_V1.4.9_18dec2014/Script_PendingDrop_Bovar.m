


clear all;
close all;

% construction de la courbe P/V pour une goutte pendante,
% Bo = 1
% (adimensionalisation par la longueur capilaire)


global RR dpdz alphac Weta Wxi verbosity beta discretization nbdim Vref;
global Bo tabP tabV tab1; 

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

tabV = [0;0;0];
tabP = [0;0;0];
tab1 = [0;0;0];


%%%
for Bo = 0.2:0.1:1.;
    Bo
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = 0.1;

[R,Z,P,V] = Newton_P(R,Z,P,V,.01);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV/10,10);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.005,100);

    if(Bo==.2)
        tabV = [tabV,[0.2;10;0.995]];
            end
        
end



figure(72);
if(size(tabV,2)>1)
plot(tabV(1,2:size(tabV,2)),tabV(3,2:size(tabV,2)),'k--');
end
hold on;
if(size(tab1,2)>1)
plot(tab1(1,2:size(tab1,2)),tab1(3,2:size(tab1,2)),'r');
end
if(size(tabP,2)>1)
plot(tabP(1,2:size(tabP,2)),tabP(3,2:size(tabP,2)),'b:');
end

figure(71);
hold on;
if(size(tabP,2)>1)
plot(tabP(1,2:size(tabP,2)),tabP(2,2:size(tabP,2)),'r');
end



tabV = [0;0;0];
tabP = [0;0;0];
tab1 = [0;0;0];


%%%
for Bo = 1:0.1:2.4;
    Bo
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = 0.3;

[R,Z,P,V] = Newton_V(R,Z,P,V,.01);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.01,100);


end

tabP = [tabP,[2.405;0;0]];

figure(72);
if(size(tabV,2)>1)
plot(tabV(1,2:size(tabV,2)),tabV(3,2:size(tabV,2)),'k--');
end
hold on;
if(size(tab1,2)>1)
plot(tab1(1,2:size(tab1,2)),tab1(3,2:size(tab1,2)),'r');
end
if(size(tabP,2)>1)
plot(tabP(1,2:size(tabP,2)),tabP(3,2:size(tabP,2)),'b:');
end

figure(71);
hold on;
if(size(tabP,2)>1)
plot(tabP(1,2:size(tabP,2)),tabP(2,2:size(tabP,2)),'r');
end

pause;

tabV = [0;0;0];
tabP = [0;0;0];
tab1 = [0;0;0];

for Bo = 2.4:.1:3.3;
    Bo
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = 0.3;

[R,Z,P,V] = Newton_V(R,Z,P,V,.01);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.01,100);


end


figure(72);
if(size(tabV,2)>1)
plot(tabV(1,2:size(tabV,2)),tabV(3,2:size(tabV,2)),'k--');
end
hold on;
if(size(tab1,2)>1)
plot(tab1(1,2:size(tab1,2)),tab1(3,2:size(tab1,2)),'r');
end
if(size(tabP,2)>1)
plot(tabP(1,2:size(tabP,2)),tabP(3,2:size(tabP,2)),'b:');
end

figure(71);
hold on;
if(size(tabP,2)>1)
plot(tabP(1,2:size(tabP,2)),tabP(2,2:size(tabP,2)),'r');
end

tabV = [0;0;0];
tabP = [0;0;0];
tab1 = [0;0;0];

for Bo = 3.3:.53/10:3.8301;
    Bo
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = 0.3;

[R,Z,P,V] = Newton_V(R,Z,P,V,.01);
[R,Z,P,V] = Loop_V(R,Z,P,V,dV,100);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.01,100);


end

tab1 = [tab1,[3.8317;0;0]];

figure(72);
if(size(tabV,2)>1)
plot(tabV(1,2:size(tabV,2)),tabV(3,2:size(tabV,2)),'k--');
end
hold on;
if(size(tab1,2)>1)
plot(tab1(1,2:size(tab1,2)),tab1(3,2:size(tab1,2)),'r');
end
if(size(tabP,2)>1)
plot(tabP(1,2:size(tabP,2)),tabP(3,2:size(tabP,2)),'b:');
end

figure(71);
hold on;
if(size(tabP,2)>1)
plot(tabP(1,2:size(tabP,2)),tabP(2,2:size(tabP,2)),'r');
end

for Bo = 3.8317:(5.1-3.8317)/12:5.1;
    Bo
R0 = Bo*t;
Z0 = 0*t;
R = R0;
Z = Z0;
P = .01;
V = 0;
dV = 0.1;

[R,Z,P,V] = Newton_V(R,Z,P,V,.01);
[R,Z,P,V] = Loop_P(R,Z,P,V,-.005,100);


end

tabV = [tabV,[5.136;0;0]];

figure(72);
if(size(tabV,2)>1)
plot(tabV(1,2:size(tabV,2)),tabV(3,2:size(tabV,2)),'k--');
end
hold on;
if(size(tab1,2)>1)
plot(tab1(1,2:size(tab1,2)),tab1(3,2:size(tab1,2)),'r');
end
if(size(tabP,2)>1)
plot(tabP(1,2:size(tabP,2)),tabP(3,2:size(tabP,2)),'b:');
end

figure(71);
hold on;
if(size(tabP,2)>1)
plot(tabP(1,2:size(tabP,2)),tabP(2,2:size(tabP,2)),'r');
end









figure(71);
xlabel('Bo');
ylabel('P');
box on;

figure(72);
xlabel('Bo');
ylabel('V');
box on;


