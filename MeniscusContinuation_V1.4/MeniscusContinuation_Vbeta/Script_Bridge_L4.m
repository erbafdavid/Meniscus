% construction de la courbe P/V pour un pont
% (adimensionalisation par le rayon)

close all;


% CAS L=4, SANS GRAVITE
global Bo RR dpdz alphac Weta Wxi verbosity beta discretization nbdim Vref;



L = 4;

nbdim = 3; % 3d or 2d
Bo = 1;
RR = 1;
dpdz = .0;
Weta = 1;
Wxi = 1;

beta = -1;
discretization = 'FD';
verbosity = 0;
Vref = pi*L;



% solution limite theorique : deux portions de sphere de rayon Rc tel que
% Rc^2 = R^2 + (Rc-L/2)^2 
F = @(x)(x^2-1-(x-L/2)^2);
Rc = fzero(F,1.5)

%alphaCC = atan(L/(4*RR))
%Rc = RR/sin(alphaCC)


%  Rc = 
Kl = 2/Rc;
Vl = 2*pi/3*(L/2)^2*(3*Rc-L/2);% formule volume d'une calotte spherique
figure(21);
plot(Kl,Vl/Vref,'go','MarkerSize',10);






% geometrie initiale : pont cylindrique
t = [0:.005:1];
R0 = RR*ones(size(t));
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
plot(P,V/Vref,'ro','MarkerSize',10);

Pbif = 1.475;

% continuation en montant P 
[R,Z,P,V] = Loop_P(R,Z,P,V,(Pbif-1)/10,10);

% point de bifurcation
[Rbif,Zbif,Pbif,Vbif] = Newton_P(R,Z,P,V,P);

[R,Z,P,V] = Loop_P(R,Z,P,V,.05,30);

%remontee
 [R,Z,P,V] = Loop_V(R,Z,P,V,.4,50);

%solution limite deux spheres qui se touchent
[Rlim,Zlim,Plim,Vlim] = Newton_P(R,Z,P,V,P);
figure(10);
plot(Rlim,Zlim,'g')




% branche asymetrique : point limite
alpha = asin(2*RR/L);
tt = pi*t-alpha;
Ra = L/2*abs(sin(tt));
Za = L/2*(((tt>0).*cos(tt)+(tt<0).*(2-cos(tt)))+cos(alpha));
Va = L^3*pi/6;
Pa = 4/L;

figure(21);
hold on;
plot(Pa,Va/Vref,'bo','MarkerSize',10);

figure(10);
plot(Ra,Za,'b');
hold on;


% construction de la branche asymetrique

%Ra = L/2*.8*abs(sin(tt))+.2*R;
%Za = L/2*(((1+tanh(5*tt))/2.*cos(tt)+(1-tanh(5*tt))/2.*(2-cos(tt)))+cos(alpha));
%Z00 = (sqrt((L/2)^2-RR^2)+L/2);
%Za = L*t;
%Ra = 1./(1+9*exp(-100*(Za-Z00).^2))+L*t.*(1-t);
%Ra = Ra-(Ra(1)-1)+t*(Ra(1)-Ra(length(Ra))) ;
%figure(10);
%plot(Ra,Za,'g:');
%hold on;
%[R,Z,P,V] = Loop_P(Ra,Za,Pa+.3,Va,+.05,100);

dpdz = -0.001;
Loop_P(Rbif,Zbif,Pbif,Vbif,(Pa-Pbif)/50,50);

dpdz = 0;


% au passage on attrape la solution d'equilibre final 
% (volume = celui des deux portions de spheres)
[Rf,Zf,Pf,Vf] = Newton_V(R0,Z0,1,0,Vl);
figure(21);
plot(Pf,Vf/Vref,'go','MarkerSize',10);
figure(10);
plot(Rf,Zf,'g')

% continuation en descendant P  pour les ponts bombes
[R,Z,P,V] = Loop_V(R0,Z0,1,Vref,.1*Vref,100);
[R,Z,P,V] = Loop_V(R,Z,P,V,.02*Vref,100);


%%%%% BRANCHE A DEUX COLS
% Resultats theoriques points de terminaison :
Rt1 = 1;
Vt1 = 2*4*pi/3;
Pt1 = 2;

%f = @(x)(4*x-2*sqrt(x^2-1))-L
Rt2 = 1.6667;
Pt2 = 2/Rt2;

% Remarque : la branche apparait pour L = 3.464 ; R = 1.15





% mise en forme des figures

figure(21);
grid on;
xlabel('P')
ylabel('V')

figure(10);
grid on;
xlabel('R')
ylabel('Z')


figure(31);
grid on;
figure(32);
grid on;

