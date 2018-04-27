% Constuction of equilibrium shapes for hanging drops, pseudo-arclength
% continuation
% pinned condition for contact line
% for Bo = 1 
%
% Script reconstructed on october 2017 to explain the usage of meniscus.m
%
% NB here m.rhog=1 and m.gamma=1 ; as a consequence the length scale is the
% capillary length. "Rhole" is the radius of the drop at the attachment location (radius of the pipe to which it is attached).
% Because of the nondimensionalization choice Rhole directly corresponds to
% the Bond number.


clear all; close all;
Rhole = 1;
m = meniscus('flatinv',100,Rhole); % initialize the meniscus with a "flat" shape
m.whichfigures = [m.whichfigures 121]% to plot thetas as function of V in figure 122
m.resetfigs; % set the figures axes and legends
m.istab = 'yes'; % to detect the bifurcation points and indicate them on fig. 10
m.discretization = 'FD'; % Finite differences (default) ; Finite elements currently works only if m.istab='no'

m = m.loop('dS',.01,1000);


