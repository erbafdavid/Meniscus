function [eta,P1,RES,lambdaminV,lambdamaxV,lambdaminP,lambdamaxP] = calceta_V(R,Z,P,alpha,ds,s0,Ka,Kb,dpdz,V)
% computation of eta : normal displacement needed to correct the curvature
% Equation : -K1(eta)+N0z*eta*dpdz = (K0 - (P+dpdz*z)) 
% Resolution as a matricial equation : A * eta = RHS

global beta discretization;

if(discretization=='FD')
    [eta,P1,RES,lambdaminV,lambdamaxV,lambdaminP,lambdamaxP] = calceta_V_FD(R,Z,P,alpha,ds,s0,Ka,Kb,dpdz,V);
else
    [eta,P1,RES,lambdaminV,lambdamaxV,lambdaminP,lambdamaxP] = calceta_V_FE(R,Z,P,alpha,ds,s0,Ka,Kb,dpdz,V);
end

