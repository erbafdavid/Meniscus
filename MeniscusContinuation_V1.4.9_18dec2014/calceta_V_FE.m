function [eta,P1,RES,lambdaminV,lambdamaxV,lambdaminP,lambdamaxP,lambdamin1,lambdamax1] = calceta_V(R,Z,P,alpha,ds,s0,Ka,Kb,dpdz,V)
% computation of eta : normal displacement needed to correct the curvature
% Equation : -K1(eta)+N0z*eta*dpdz = (K0 - (P+dpdz*z)) 
% Resolution as a matricial equation : A * eta = RHS
global beta verbosity vecPGG;

V = V/pi;

A = 0;
N = length(R);

for i=2:N-1
%  - (d etabar / d s0 * d eta /d s0) *R
A(i,i-1) = (R(i)+R(i-1))/2*1/ds(i-1);
A(i,i) = - R(i)*(1/ds(i)+1/ds(i-1));
%if (i~=N-1) 
A(i,i+1) = (R(i)+R(i+1))/2*1/ds(i); 
%end
%  (Ka^2 + Kb^2 + dpdz*cos(alpha)) etabar eta 
A(i,i-1) = A(i,i-1)  + (  (Ka(i)^2+Kb(i)^2+dpdz*cos(alpha(i)))*R(i) ...
                        + (Ka(i-1)^2+Kb(i-1)^2+dpdz*cos(alpha(i-1)))*R(i-1) )/2 * (ds(i-1)/6);  
A(i,i) = A(i,i) + (Ka(i)^2+Kb(i)^2+ dpdz*cos(alpha(i))) * (R(i)*(ds(i-1)+ds(i))/3+(ds(i)^2+ds(i-1)^2)*cos(alpha(i))/12 );
%if (i~=N-1) 
A(i,i+1) = A(i,i+1)  + (  (Ka(i)^2+Kb(i)^2+dpdz*cos(alpha(i)))*R(i) ...
                        + (Ka(i+1)^2+Kb(i+1)^2+dpdz*cos(alpha(i+1)))*R(i+1) )/2 * (ds(i)/6); 
%end

end

% Right-hand side
for i=2:N-1
    RHS(i) = (Ka(i)+Kb(i)-P+dpdz*Z(i))*R(i)*(ds(i-1)+ds(i))/2;
end
   

% case i=1
if (abs(R(1)) < 1e-6) %for bubble/drop : point 1 is on symmetry axis
  A(1,1) = -1/ds(1)*(R(2)/2); 
  A(1,2) = A(2,1);
  A(1,1) = A(1,1)+    (Ka(1)^2+Kb(1)^2+dpdz*cos(alpha(1))) * (ds(1)^2/12);
  RHS(1) = (Ka(1)+Kb(1)-P+dpdz*Z(1))*ds(1)^2/6; %%% int ( etatilde(1) * (K-P) R ds)
else % for bridge : point 1 is pinned
  A(1,1) = 1;
  A(1,2) = 0;
  RHS(1) = 0;  
end


if (beta<0) 
    % "pinned" case
A(N,N) = -1e12;
A(N,N-1) = 0;
A(N-1,N)=0;
RHS(N) = 0;
else % fixed-angle case
  A(N,N) = -1/ds(N-1)*R(N); 
  A(N,N-1) = A(N-1,N);
  A(N,N) = A(N,N) + (   (Ka(N)^2+Kb(N)^2+dpdz*cos(alpha(N)))*R(N) ...
                      + (Ka(N-1)^2+Kb(N-1)^2+dpdz*cos(alpha(N-1)))*R(N-1) )/2 * (ds(N-1)/3); 
  RHS(N) = -alpha(N)+beta;
end


% last equation : volume variation
for i=2:N-1
A(N+1,i) =  R(i)*(ds(i)+ds(i-1))/2; % test sign
end;
if (R(1)~=0)
A(N+1,1) = R(2)/2*ds(1)/2; % test sign
end
if(beta<0)
   A(N+1,N) = R(N-1)/2*ds(N-1)/2; 
end
A(N+1,N+1) = 0;

% last column : effet of dP
for i=1:N-1
A(i,N+1) = A(N+1,i);
end


% Last equation : volume
Vcurrent = 0;
 for i=2:N
 Vcurrent = Vcurrent + ( (R(i)+R(i-1))^2/4 + (R(i)-R(i-1))^2/12 )*(Z(i)-Z(i-1));
 end;
RHS(N+1) = (V-Vcurrent)/2; % test sign

if (verbosity >2)
Ad = diag(A(1:N-1,1:N-1))'
Ainf = diag(A(2:N,1:N-1))'
Asup = diag(A(1:N-1,2:N))'
Alast = A(N,:)
RHS
end


%Resolution

etaS = (A\RHS')';
P1 = etaS(N+1);
eta = etaS(1:N);
RES = max(abs(RHS));

 
% eigenvalues of the volume-controled problem

A(N,N-1) = A(N-1,N);
A(1,2) = A(2,1);

EE = eig(A);
lambdaminV = EE(abs(EE)==min(abs(EE)));
lambdamaxV = EE(real(EE)==max(real(EE)));
lambdaminV=lambdaminV(1);
lambdamaxV = min(lambdamaxV(1),1);



%RHS(N) = RHS(N)/TGVvol; % test

    
% eigenvalues of the pressure-controled problem

EE = eig(A(1:N,1:N));
lambdaminP = EE(abs(EE)==min(abs(EE)));
lambdamaxP = EE(real(EE)==max(real(EE)));
lambdaminP = lambdaminP(1);
lambdamaxP = min(lambdamaxP(1),1);



if (abs((V*pi-1.399)*(V*pi-1.6603))<1e-4) %% DEBUG FOR FE/FD comparison
[Vec,D] = eig(A(1:N,1:N));
VecP = Vec(:,abs(EE)==min(abs(EE)));

lambdaminP
vecPGG
vecPGG(N)
RESFE = (A(1:N,1:N)*vecPGG)'
RESFE(N)

figure(67);
hold on;
plot(s0,VecP);
figure(68);hold on;
plot(s0,VecP,'r',s0,-vecPGG,'b');
figure(69);hold on;
plot(s0,RESFE,'r');

end;


% eigenvalues of the m=1 problem

AA1 = A(1:N,1:N);
for i=2:N-1
if(R(i-1)~=0)
    AA1(i,i-1) = AA1(i,i-1)  + (  1/R(i)  + 1/R(i-1) )/2 * (ds(i-1)/6); 
else
    AA1(i,i-1) = AA1(i,i-1)  + (  2/R(i)  ) * (ds(i-1)/6);  % valeur de AA1(2,1) sur l'axe
end   
AA1(i,i) = AA1(i,i) + 1/R(i)*(ds(i-1)+ds(i))/3; 
%if (i~=N-1) 
AA1(i,i+1) = AA1(i,i+1)  + ( 1/R(i) +1/R(i+1) )/2 * (ds(i)/6); 
end    

AA1(1,1) = 1; AA1(1,2) = 0;
    
if(beta<0)
AA1(N,N) = 1; AA1(N,N-1) = 0;
else
AA1(N,N) = AA1(N,N) + 1/R(N)*(ds(N))/3; 
AA1(N,N-1) = AA1(N-1,N);
end


EE = eig(AA1);
lambdamin1 = EE(abs(EE)==min(abs(EE)));
lambdamax1 = EE(real(EE)==max(real(EE)));
lambdamin1 = lambdamin1(1);
lambdamax1 = min(lambdamax1(1),1);





