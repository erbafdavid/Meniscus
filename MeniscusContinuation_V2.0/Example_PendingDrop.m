% Constuction of equilibrium shapes for hanging drops, pseudo-arclength
% continuation
% pinned condition for contact line
% for Bo = 1 
%
% Script reconstructed on october 2017 to explain the usage of meniscus.m
%
% nb here m.rhog=1 and m.gamma=1 so the Bond number is Rhole


clear all; close all;
Rhole = 1;
m = meniscus('flat',100,Rhole); % initialize the meniscus
resetfigs(m); % set the figures axes and legends
m = loop(m,'S',.01,1000);


