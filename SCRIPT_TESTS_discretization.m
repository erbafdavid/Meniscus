clear all; close all;

% tests de la forme d'equilibre dans le cas 'pined'

m=meniscus('sphere',50,[.5,pi/2]); % approximate shape = semisphere
m.resetfigs; 
m.rhog = -1; % because rhog is actually (rho1-rho2) g and rho2 > rho1
m.typeend = 'pined';
m.idebug = 10;
m.discretization = 'FD';
m = step(m,'V',-0.1); % compute the shape with fixed volume
m = step(m,'P',-0.01); % compute the shape with fixed pressure





pause;

% test of Finite elements method
m=meniscus('sphere',50,[.5,pi/2]); % approximate shape = semisphere
m.resetfigs; 
m.rhog = -1; % because rhog is actually (rho1-rho2) g and rho2 > rho1
m.typeend = 'pined';
m.idebug = 10;
m.discretization = 'FE';

%test of step p with convenient starting point for P

m.P = 4.33;
m = step(m,'P',0.000); % re compute the shape

% test successful

%test of step V with convenient starting point for P

m = step(m,'V',0.01); % re compute the shape

% test successful !!!!


% next step : angle imposé !

