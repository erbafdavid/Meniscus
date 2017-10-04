function [eta,RES,lambdamin,lambdamax] = calceta_P(R,Z,alpha,ds,s0,Ka,Kb,dpdz,P)
% computation of eta : normal displacement needed to correct the curvature
% Equation : -K1(eta)+N0z*eta = (K0 - (P-z)) 
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
A(i,i) = A(i,i) - dpdz*cos(alpha(i));
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
  A(1,1) = A(1,1)+(Ka(1)^2+Kb(1)^2)-dpdz;
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

A(2,2)

if (A(2,2)*0==0) % to check if NaN

EE = eig(A);
lambdaminP = EE(abs(EE)==min(abs(EE)));
lambdamaxP = EE(real(EE)==max(real(EE)));
lambdaminP = lambdaminP(1);
lambdamaxP = lambdamaxP(1);

if (verbosity >2)
    A
Ad = diag(A)
Ainf = diag(A(2:N,1:N-1))
Asup = diag(A(1:N-1,2:N))
Alast = A(N,:)
RHS
end

eta = (A\RHS')';

    
RES = max(abs(RHS));


lambdamin = lambdaminP;
lambdamax = lambdamaxP;

else
    
eta = 0.*R;
RES = 5000;
lambdamin = 1;
lambdamax = 1;

end


