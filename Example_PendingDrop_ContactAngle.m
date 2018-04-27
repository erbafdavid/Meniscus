% Constuction of equilibrium shapes for hanging drops,  
% with fixed contact angle
%
% Script reconstructed on october 2017 to explain the usage of meniscus.m
%
% nb here m.rhog=1 and m.gamma=1 so the length scale is the capillary
% length


clear all; close all;
beta = 60; % angle in degrees
Rinit = .4; 
m = meniscus('capinv',200,[Rinit,beta]); % initialize the meniscus
m.whichfigures = [m.whichfigures 122];
m.resetfigs; % set the figures axes and legends
%m.istab = 'yes';
m.discretization = 'FD'; % Finite differences (default) ; Finite elements currently works only if m.istab='no'

% first arclencght-continuation loop with coarse step
m = m.loop('dS',.05,200);




% after 85 steps the Newton fails 

%so we continue with a smaller step to take the turn of the curve
%m = m.loop('dS',.005,50);

%after the turn is successfully taken we continue with a larger step
%m = m.loop('dS',.05,100);%


