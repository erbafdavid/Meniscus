
%% Construction of a family of shapes for
%% HANGING DROPS with PINNED CONTACT LINES

close all;

thetas = 30; % degrees;
m = meniscus('cap',200,[.5,thetas]);
m.discretization = 'FD' ; % finite differences
m.istab = 'yes';
m.Vref = 0.3;% marche mieux pour les petits angles 

m.whichfigures = [10 20 28];
m.resetfigs;


m = m.loop('dS',.05,50);
m = m.loop('dS',.1,50);
m = m.loop('dS',1,50);
m = m.loop('dS',1,50);
m = m.loop('dS',1,50);
m = m.loop('dS',1,50);