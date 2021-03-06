% Constuction of equilibrium shapes for sessile drops,  
% with pinned contact line
%
% Script reconstructed on october 2017 to explain the usage of meniscus.m
%
% This case uses DIMENSIONAL for a water droplet on PMMS with contact angle
% 85 degrees (average static angle) ; 78.5 degrees (recessing angle) ; 
% 92 degrees (advancing angle).
%
% NOTE : we use cm/g/s units instead of international units (m/kg/s)
% otherwise the matrices are badly scalled.
% 
% gamma = 0.0728 N/m = 72.8 g/s^2
% rho g = 1000 kg/m^3 * 9.81 m/s^2 = 981 g/cm^2/s^2


clear all; close all;
Rinit = 0.3; 
beta = 85; 
m = meniscus('cap',200,[Rinit,beta]); % initialize the meniscus
m.plotmeniscus(10,'k:')
m.gamma = 72.8; 
m.rhog = 981;
m.typestart = 'angle';
m.discretization = 'FD'; % Finite differences (default) ; Finite elements currently works only if m.istab='no'
m.whichfigures = [m.whichfigures 122];
m.resetfigs; % set the figures axes and legends
m = m.step('dV',0); % first computation

m = m.loop('dV',0.05,100);


beta = 78.5; 
m = meniscus('cap',200,[Rinit,beta]); % initialize the meniscus
m.plotmeniscus(10,'k:')
m.gamma = 72.8; 
m.rhog = 981;
m.typestart = 'angle';
m.discretization = 'FD'; % Finite differences (default) ; Finite elements currently works only if m.istab='no'
m.whichfigures = [m.whichfigures 122];
m = m.step('dV',0); % firt computation

%m = m.loop('dV',0.05,20);
%m = m.loop('dV',0.05,20)
%m = m.loop('dV',0.05,20)
%m = m.loop('dV',0.05,20)
%m = m.loop('dV',0.05,20)


m = m.loop('dV',0.05,100);




beta = 92; 
m = meniscus('cap',200,[Rinit,beta]); % initialize the meniscus
m.plotmeniscus(10,'k:')
m.gamma = 72.8; 
m.rhog = 981;
m.typestart = 'angle';
m.discretization = 'FD'; % Finite differences (default) ; Finite elements currently works only if m.istab='no'
m.whichfigures = [m.whichfigures 122];
m = m.step('dV',0); % firt computation

%first arclencght-continuation loop with coarse step
%m = m.loop('dS',.05,200);

m = m.loop('dV',0.05,100);



