function [eta,RES,lambdamin,lambdamax] = calceta_P(R,Z,alpha,ds,s0,Ka,Kb,dpdz,P)
% computation of eta : normal displacement needed to correct the curvature
% Equation : -K1(eta)+N0z*eta = (K0 - (P-z)) 

global beta discretization;

if(discretization=='FD')
    [eta,RES,lambdamin,lambdamax] = calceta_P_FD(R,Z,alpha,ds,s0,Ka,Kb,dpdz,P);
else
    [eta,RES,lambdamin,lambdamax] = calceta_P_FE(R,Z,alpha,ds,s0,Ka,Kb,dpdz,P);
end

