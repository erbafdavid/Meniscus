% script pour calculer le menisque du cas "Viola, Brun & Gallaire"
% cas Bo = 500

m = meniscus('flat',200,1);
m.gamma = 1/500;
m.discretization = 'FD'; % case FE with gamma different from 1 still to be tested
m = m.step('P',-0.035),plotmeniscus(m,10,'k:');

anglecontact = m.alpha(1)*180/pi
