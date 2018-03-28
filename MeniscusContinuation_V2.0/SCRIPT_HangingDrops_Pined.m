
%% Construction of a family of shapes for
%% HANGING DROPS with PINNED CONTACT LINES

close all;
m = meniscus('flatinv',100,1);
m.istab = 'yes'; % to detect the bifurcation points and indicate them on fig. 10

m = m.loop('dS',.02,200);
