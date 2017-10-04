function [eta,P1,RES,lambdaminV,lambdamaxV,lambdaminP,lambdamaxP] = calceta_V(R,Z,P,alpha,ds,s0,Ka,Kb,dpdz,V)
% computation of eta : normal displacement needed to correct the curvature
% Equation : -K1(eta)+N0z*eta*dpdz = (K0 - (P+dpdz*z)) 
% Resolution as a matricial equation : A * eta = RHS
global beta discretization;

V = V/pi;

A = 0;
N = length(R);

for i=2:N-1
% (d etabar / d s0 * d eta /d s0) *R
A(i,i-1) = (R(i)+R(i-1))/2*1/ds(i-1);
A(i,i) = - R(i)*(1/ds(i)+1/ds(i-1));
if (i~=N-1) 
A(i,i+1) = (R(i)+R(i+1))/2*1/ds(i); 
end
%  (Ka^2 + Kb^2- dpdz*cos(alpha)) etabar eta 
A(i,i-1) = A(i,i-1)  + (  (Ka(i)^2+Kb(i)^2+dpdz*cos(alpha(i)))*R(i) ...
                        + (Ka(i-1)^2+Kb(i-1)^2+dpdz*cos(alpha(i-1)))*R(i-1) )/2 * (ds(i-1)/6);  
A(i,i) = A(i,i) + (Ka(i)^2+Kb(i)^2+ dpdz*cos(alpha(i)))*R(i)*(ds(i-1)+ds(i))/3; 
%if (i~=N-1) 
A(i,i+1) = A(i,i+1)  + (  (Ka(i)^2+Kb(i)^2+dpdz*cos(alpha(i)))*R(i) ...
                        + (Ka(i+1)^2+Kb(i+1)^2+dpdz*cos(alpha(i+1)))*R(i+1) )/2 * (ds(i)/6); 
%end

end

% Right-hand side
for i=2:N-1
    RHS(i) = (Ka(i)+Kb(i)-P+dpdz*Z(i))*R(i)*(ds(i-1)+ds(i))/2;
end
    RHS(1) = 0; 

% case i=1
if (abs(R(1)) < 1e-6) %for bubble/drop : point 1 is on symmetry axis
  A(1,1) = -1/ds(1)*(R(2)/2); 
  A(1,2) = A(2,1);
  A(1,1) = A(1,1) + (   (Ka(1)^2+Kb(1)^2+dpdz*cos(alpha(1)))*R(1) ...
                      + (Ka(2)^2+Kb(2)^2+dpdz*cos(alpha(2)))*R(2) )/2 * (ds(1)/3); 
  R(1) = 0; % in case of trouble...
else % for bridge : point 1 is pinned
  A(1,1) = 1;
  A(1,2) = 0;
  RHS(1) = 0;  
end


% case i=N (pinned)
A(N,N) = 1;
RHS(N) = 0;



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


%A

EE = eig(A);
lambdaminV = EE(abs(EE)==min(abs(EE)));
lambdamaxV = EE(real(EE)==max(real(EE)));
lambdaminV = lambdaminV(1);
lambdamaxV = lambdamaxV(1);



%EE




etaS = (A\RHS')';

P1 = etaS(N+1); % test sign
eta = etaS(1:N);

%V
%Vcurrent




RES = max(abs(RHS));
    

EE = eig(A(1:N,1:N));
lambdaminP = EE(abs(EE)==min(abs(EE)));
lambdamaxP = EE(real(EE)==max(real(EE)));
lambdaminP=lambdaminP(1);
lambdamaxP=lambdamaxP(1);

V = V*pi;

%diag(A)
%diag(A(1:N-1,2:N))
%diag(A(2:N,1:N-1))

%lambdaminP 
%lambdamaxP 
%lambdaminV
%lambdamaxV


%lambdamax = lambdamaxV;



