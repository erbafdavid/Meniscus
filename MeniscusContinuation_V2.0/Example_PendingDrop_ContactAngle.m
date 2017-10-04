% Constuction of equilibrium shapes for hanging drops,  
% with fixed contact angle
%
% Script reconstructed on october 2017 to explain the usage of meniscus.m
%
% nb here m.rhog=1 and m.gamma=1 so the Bond number is Rhole


clear all; close all;
beta = pi/2;
Rinit = .4; 
m = meniscus('sphere',200,[Rinit,beta]); % initialize the meniscus
resetfigs(m); % set the figures axes and legends

% first arclencght-continuation loop with coarse step
m = loop(m,'S',.05,100);
% after 85 steps the Newton fails 

%so we continue with a smaller step to take the turn of the curve
m = loop(m,'S',.005,50);

%after the turn is successfully taken we continue with a larger step
m = loop(m,'S',.05,100);


