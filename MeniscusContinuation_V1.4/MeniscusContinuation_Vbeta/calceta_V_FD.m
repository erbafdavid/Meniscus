function [eta,P1,RES,lambdaminV,lambdamaxV,lambdaminP,lambdamaxP,lambdamin1,lambdamax1] = calceta_V(R,Z,P,alpha,ds,s0,Ka,Kb,dpdz,V)
% computation of eta : normal displacement needed to correct the curvature
% Equation : -K1(eta)+N0z*eta*dpdz = (K0 - (P+dpdz*z)) 
% Resolution as a matricial equation : A * eta = RHS
global beta verbosity nbdim;


A = 0;
N = length(R);

for i=2:N-1
% d^2 eta / d^2 s0
A(i,i-1) = 2*ds(i)/(ds(i)*ds(i-1)*(ds(i)+ds(i-1)));  
A(i,i) = -2/(ds(i)*ds(i-1));
%if (i~=N-1)
A(i,i+1) = 2*ds(i-1)/(ds(i)*ds(i-1)*(ds(i)+ds(i-1)));  
%end
% T0r/r* d eta / d s0
if (nbdim==3)
    A(i,i-1) = A(i,i-1) -ds(i)^2/((ds(i)+ds(i-1))*ds(i)*ds(i-1))*cos(alpha(i))/R(i);
    A(i,i) = A(i,i) +  (ds(i)^2-ds(i-1)^2)/((ds(i)+ds(i-1))*ds(i)*ds(i-1))*cos(alpha(i))/R(i);
    %if (i~=N-1)
        A(i,i+1) = A(i,i+1) +ds(i-1)^2/((ds(i)+ds(i-1))*ds(i)*ds(i-1))*cos(alpha(i))/R(i);
    %end
end
%  (Ka^2 + Kb^2) eta
A(i,i) = A(i,i) + 1*(Ka(i)^2+Kb(i)^2); %tentative marche pas mieux
% + N0z*eta*dpdz%
A(i,i) = A(i,i) + dpdz*cos(alpha(i));
end;

% Right-hand side
for i=2:N-1
    RHS(i) = Ka(i)+Kb(i)-P+dpdz*Z(i);
end


% case i=1
if (abs(R(1)) < 1e-6) %for bubble/drop : point 1 is on symmetry axis
  if (nbdim==3)
    A(1,1) = -4/(ds(1)^2);  %%% facteur 4 car 1/R d eta / ds = d^2 eta / d s^2 
    A(1,2) = 4/(ds(1)^2);
  else
    A(1,1) = -2/(ds(1)^2);  %%% facteur 4 car 1/R d eta / ds = d^2 eta / d s^2 
    A(1,2) = 2/(ds(1)^2);
  end
  A(1,1) = A(1,1)+(Ka(1)^2+Kb(1)^2)+dpdz;
  RHS(1) = Ka(1)+Kb(1)-P+dpdz*Z(1);
  R(1) = 0; % in case of trouble...
else % for bridge : point 1 is pinned
  A(1,1) = 1;
  A(1,2) = 0;
  RHS(1) = 0;  
end


if (beta<0) 
    % "pinned" case
A(N,N) = 1;
RHS(N) = 0;
else
   A(N,N-2) =  1/(2*ds(N-1));
   A(N,N-1) = -4/(2*ds(N-1));
   A(N,N)   =  3/(2*ds(N-1));

%A(N,N-1) = -1/(ds(N-1));
%A(N,N)   =  1/(ds(N-1));

   RHS(N)   =  (alpha(N)-beta);
end

TGVvol = 1; %%% tentative d'ameliorer la convergence => marche pas tellement

% last equation : volume variation
  % WARNING works only in 3D !

for i=2:N-1
A(N+1,i) =  (ds(i-1)/2*(2*R(i)+R(i-1))/3 + ds(i)/2*(2*R(i)+R(i+1))/3)*2*pi*TGVvol;
end;
%if (R(1)~=0)
A(N+1,1) = ds(1)/2*(2*R(1)+R(2))/3*2*pi*TGVvol;
%end
A(N+1,N) = ds(N-1)/2*(2*R(N)+R(N-1))/3*2*pi*TGVvol;

% last column : effet of dP
for i=2:N-1
A(i,N+1) = +1;
end;
if(R(1)==0) 
    A(1,N+1) = +1;
end

%if(beta>0) 
%    A(N,N+1) = +1;
%end

% Last equation : volume
Vcurrent = 0;
 for i=2:N
 Vcurrent = Vcurrent + ( (R(i)+R(i-1))^2/4 + (R(i)-R(i-1))^2/12 )*(Z(i)-Z(i-1))*pi;
 end;
RHS(N+1) = (V-Vcurrent)*TGVvol;


%Resolution

etaS = (A\RHS')';
P1 = etaS(N+1);
eta = etaS(1:N);
RES = max(abs(RHS));

% eigenvalues of the volume-controled problem

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

% eigenvalues of the m=1 problem

AA1 = A(1:N,1:N);
for i=2:N
    AA1(i,i) = AA1(i,i) - 1/R(i)^2;
end
AA1(1,1) = 1; AA1(1,2) = 0;
AA1(N,N) = 1; AA1(N,N-1) = 0;



EE = eig(AA1);
lambdamin1 = EE(abs(EE)==min(abs(EE)));
lambdamax1 = EE(real(EE)==max(real(EE)));
lambdamin1 = lambdamin1(1);
lambdamax1 = min(lambdamax1(1),1);




