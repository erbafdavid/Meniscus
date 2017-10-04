function [eta,RES,lambdamin,lambdamax] = calceta_P(R,Z,alpha,ds,s0,Ka,Kb,dpdz,P)
% computation of eta : normal displacement needed to correct the curvature
% Equation : -K1(eta)+N0z*eta = (K0 - (P-z)) 
% Resolution as a matricial equation : A * eta = RHS
global beta;
% on peut virer Bo dans les arguments

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
A(i,i-1) = A(i,i-1)  + (  (Ka(i)^2+Kb(i)^2-dpdz*cos(alpha(i)))*R(i) ...
                        + (Ka(i-1)^2+Kb(i-1)^2-dpdz*cos(alpha(i-1)))*R(i-1) )/2 * (ds(i-1)/6);  
A(i,i) = A(i,i) + (Ka(i)^2+Kb(i)^2- dpdz*cos(alpha(i)))*R(i)*(ds(i-1)+ds(i))/3; 
if (i~=N-1) 
A(i,i+1) = A(i,i+1)  + (  (Ka(i)^2+Kb(i)^2-dpdz*cos(alpha(i)))*R(i) ...
                        + (Ka(i+1)^2+Kb(i+1)^2-dpdz*cos(alpha(i+1)))*R(i+1) )/2 * (ds(i)/6); 
end

end



% Right-hand side
for i=2:N-1
    RHS(i) = (Ka(i)+Kb(i)-P+dpdz*Z(i))*R(i)*(ds(i-1)+ds(i))/2;
end
    RHS(1) = 0; 

% case i=1
if (abs(R(1)) < 1e-6) %for bubble/drop : point 1 is on symmetry axis
  A(1,1) = -1/ds(1)*(R(2)/2);  %%% facteur 4 car 1/R d eta / ds = d^2 eta / d s^2 (mais en fait 2...) 
  A(1,2) = A(2,1);
  A(1,1) = A(1,1)+ (   (Ka(1)^2+Kb(1)^2-dpdz*cos(alpha(1)))*R(1) ...
                     + (Ka(2)^2+Kb(2)^2-dpdz*cos(alpha(2)))*R(2) )/2 * (ds(1)/3); 
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
end
    


EE = eig(A);
lambdaminP = EE(abs(EE)==min(abs(EE)));
lambdamaxP = EE(real(EE)==max(real(EE)));
lambdaminP =lambdaminP(1);
lambdamaxP = lambdamaxP(1);



eta = (A\RHS')';

%eta(N) = 0;
    
RES = max(abs(RHS));



lambdamin = lambdaminP;
lambdamax = lambdamaxP;

