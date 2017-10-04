function [alpha,ds,s0,Ka,Kb] = calcgeom(R,Z)
% calculation of geometrical properties of a curve in the (R,Z) plane
% input : (R,Z) two arrays of points describing the curve (length N)
%   starting from the nose of the drop (R(1) = 0)
%   ending on the nozzle (R(N) = BOND NUMBER)
% output : 5 vectors of length N-1 (except s0 and alpha : N)
% alpha : angle of tangent 
% ds, s0 : curvilinear abcisse (step,total)
% Ka, Kb : pricipal curvature radii

global nbdim;

N = length(R);


%%%% FIRST METHOD

ds(1) = sqrt((R(2)-R(1))^2+(Z(2)-Z(1))^2);
s0(1) = 0;
%alpha(1) = 0;
%Ka(1) = curv3P(-R(2),Z(2),R(1),Z(1),R(2),Z(2));
%Kb(1) = Ka(1);

for i=2:N-1
    ds(i) = sqrt( (R(i+1)-R(i))^2+(Z(i+1)-Z(i))^2);
    s0(i) = s0(i-1)+ds(i-1);
%    Ka(i) = curv3P(R(i-1),Z(i-1),R(i),Z(i),R(i+1),Z(i+1));
%    alpha(i) = atan2(Z(i)-Z(i-1),R(i)-R(i-1))+asin(ds(i)*Ka(i)/2);
%    Kb(i) = sin(alpha(i))/R(i);
end

     s0(N) = s0(N-1)+ds(N-1);
     
%    alpha(N) = atan2(Z(N)-Z(N-1),R(N)-R(N-1))+asin(ds(N-1)*Ka(N-1)/2);
%    Ka(N) = Ka(N-1);%%% curv3P(R(N-1),Z(N-1),R(N),Z(N),1,0);
%    Kb(N) = sin(alpha(N))/R(N);

    
    %%%% SECOND METHOD
    for i=2:N-1
        alpha(i) =  (ds(i)*atan2(Z(i)-Z(i-1),R(i)-R(i-1)) +  ds(i-1)*atan2(Z(i+1)-Z(i),R(i+1)-R(i)))/(ds(i-1)+ds(i));
        Ka1(i) = (atan2(Z(i+1)-Z(i),R(i+1)-R(i))-atan2(Z(i)-Z(i-1),R(i)-R(i-1)))/((ds(i)+ds(i-1))/2);
        if (nbdim == 3) 
            Kb1(i) = sin(alpha(i))/R(i);
        else
            Kb1(i) = 0;
        end
    end;
    
    % i = N : pinned point (geom. quantities only needed for visualization and for angle-fixed case)
    
  %  alpha(N) =2*alpha(N-1)-alpha(N-2)
     alpha(N) = alpha(N-2)+2*ds(N-1)*(Ka1(N-1)+Kb1(N-1)-sin(alpha(N-1))/R(N-1));
    Ka1(N) =2*Ka1(N-1)-Ka1(N-2); 
    Kb1(N) =2*Kb1(N-1)-Kb1(N-2); 
    
  %  alpha(N-2:N)
    
    if(abs(R(1))>1e-6)
    %  i = 1 : pinned point
    alpha(1) = 2*alpha(2)-alpha(3); 
    Ka1(1) = 2*Ka1(2)-Ka1(3); 
    Kb1(1) = 2*Kb1(2)-Kb1(3); 

    else
     %  i = 1 : axis   
     alpha(1) = 0.;
    Ka1(1) = 2*atan2(Z(2)-Z(1),R(2)-R(1) )/ds(1); % facteur 2 rajouté 
      if (nbdim == 3) 
          Kb1(1) = Ka1(1);
      end
    end     
    
    %%%% DEBUG
    verb = 1;
    if(verb >10)
    figure(430);
    plot(s0,alpha,'r+',s0,alpha1,'b+')
    figure(431);
    plot(s0,Ka,'r+',s0,Ka1,'b+')
     figure(432);
    plot(s0,Kb,'r+',s0,Kb1,'b+')
    figure(433);
    plot(s0,Ka+Kb+Z,'r+',s0,Ka1+Kb1+Z,'b-')
    figure(434);
    plot(s0,Ka1-Ka,'r+',s0,Kb1-Kb,'b-')
    end
    
     if(verb >8)
   figure(435);
   hold on;
   ds(N) = ds(N-1);
    plot(s0,ds,'r-+')
    end
    
    % test
    Ka = Ka1;
    Kb = Kb1;
        % Ca n'arrange pas le probleme !
    