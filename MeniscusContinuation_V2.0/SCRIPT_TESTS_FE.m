
%% SCRIPT TO TEST THE OPERATION OF 'MENISCUS.m'
%% USING FINITE DIFFERENCES.

%% First series of tests : construction of a family of shapes for
%% SESSILE DROPS with PINNED CONTACT LINES
%% This is initiated by the a 'flat' initial condition (numbered in the
%% direct sense, namely point 1 is at the pinning point and point N on the axis

close all;
m = meniscus('flat',100,1);
m.rhog=1;% normally this is default value
m.idebug = 0;% set to 100 for instance to enter "debug mode"
m.discretization = 'FE'; % idem;

% basis tests in 'P' mode and 'V' mode
m.step('P',0.1);plotmeniscus(m,10,'b--');pause;
m.step('V',0.1);plotmeniscus(m,10,'r--');pause;

% we now scan the family of shapes using 'dS', 'dP' and 'dV' modes to test the
% operation of each
m = m.loop('dP',0.05,10);
m = m.loop('dV',0.05,10);
%m = m.loop('dS',.05,50);
pause;


%% SECOND series of tests : construction of a family of shapes for
%% HANGING DROPS with PINNED CONTACT LINES
%% This is initiated by the a 'flatinv' initial condition (numbered in the
%% reveresed sense, namely point 1 is at the axis and point N at the pinning point

close all;
m = meniscus('flatinv',100,1);
m.rhog=1;% normally this is default value
m.idebug = 0;% set to 100 for instance to enter "debug mode"
m.discretization = 'FE'; % idem;

% basis tests in 'P' mode and 'V' mode
m.step('P',0.1);plotmeniscus(m,10,'b--');pause;
m.step('V',0.1);plotmeniscus(m,10,'r--');pause;

% we construct the shape using 'dS', 'dP' and 'dV' modes to test the
% operation of each
m = m.loop('dP',0.05,10);
m = m.loop('dV',0.05,10);
%m = m.loop('dS',.02,200);

%% THIRD SERIES OF TESTS :
%% SESSILE DROPS WITH IMPOSED CONTACT LINE


close all;
R0 = 0.2; % we start with a small droplet, so that it is almost spherical
thetas = 60; % (in degrees) 
m = meniscus('cap',100,[.2,thetas]);
m.discretization = 'FE';
% tests in P and V modes
m.step('dV',0.0);plotmeniscus(m,10,'b--');pause;
m.step('dP',0.01);plotmeniscus(m,10,'r--');pause;
m = m.loop('dP',0.05,10);
m = m.loop('dV',0.05,10);
%m = m.loop('dS',.05,50);
pause;

%% FOURTH SERIES OF TESTS :
%% HANGING DROPS WITH IMPOSED CONTACT LINE


close all;
R0 = 0.2; % we start with a small droplet, so that it is almost spherical
thetas = 60; % (in degrees) 
m = meniscus('capinv',100,[.2,thetas]);
m.discretization = 'FE';
% tests in P and V modes
m.step('dV',0.0);plotmeniscus(m,10,'b--');pause;
m.step('dP',0.01);plotmeniscus(m,10,'r--');pause;
m = m.loop('dP',0.05,10);
m = m.loop('dV',0.05,10);
%m = m.loop('dS',.05,50);
pause;


