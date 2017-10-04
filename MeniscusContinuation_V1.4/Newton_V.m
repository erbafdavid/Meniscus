function [R,Z,P,V,conv] = Newton_V(R0,Z0,P0,V0,newV)
% computation of equilibrium solution through Newton method
% input : 
% [R0,Z0] "guess" to start the iteration
% Bo Bond number ( radius divided by capilary length )
% P : nondimensional pressure at nozzle (height of liquid column divided by capilary length)
% output : 
% [R,Z] curve for equilibrium solution
% V  volume

global Bo RR dpdz alphac Weta Wxi verbosity beta discretization Vref; 

V = newV;

N = length(R0);
 
  R = R0;
  Z = Z0;
  P = P0;
  RES = 1;
  epsilon = 1e-6;
  it = 0;
  itmax = 100;

%  figure(11);
%  plot(R,Z,'r+-');
%  hold on;


  if(verbosity > 2)
    figure(24);
    plot(R,Z,'r');
    hold on;
  end
  
  while ((RES > epsilon)&&(it<itmax)&&(RES<100))
      [alpha,ds,s0,Ka,Kb] = calcgeom(R,Z);
      if(max(abs(ds))>1)
          RES = 20000;
      end
      [eta,P1,RES,lambdaminV,lambdamaxV,lambdaminP,lambdamaxP] = calceta_V(R,Z,P,alpha,ds,s0,Ka,Kb,dpdz,V);
      
%      R = R + Weta*eta.*sin(alpha);
%      Z = Z - Weta*eta.*cos(alpha);
%       if (beta<0)
%          xiN = 0;
%       else
%          [alpha,ds,s0,Ka,Kb] = calcgeom(R,Z);
%          xiN = -Z(N)/sin(alpha(N));
%      end
%      xi = calcxi(R,Z,0,xiN);
%      R = R  + Wxi*xi.*cos(alpha);
%      Z = Z  + Wxi*xi.*sin(alpha);
      
 %    figure(33);
 %    plot(R,Z,'r+-')
 %    hold on;
    
      R = R + Weta*eta.*sin(alpha);
      Z = Z - Weta*eta.*cos(alpha);
 %    figure(33);
 %     plot(R,Z,'bo-')
    
         if (beta<0)
          xiN = 0;
          ZN = 0;
       else
          [alpha,ds,s0,Ka,Kb] = calcgeom(R,Z);
%          xiN = -Z(N)/sin(alpha(N)); ZN = 0; % methode 1 : marche pour 45<beta<135
          xiN = 0;   ZN = Z(N); %methode 2
      end
      xi = calcxi(R,Z,0,xiN);
      R = R  + Wxi*xi.*cos(alpha);
      Z = Z  + Wxi*xi.*sin(alpha)-ZN*(Z-Z(1))/(Z(N)-Z(1));
 %     figure(33);
 %      plot(R,Z,'go-')
  %    pause ;

      
      P = P + P1*1;
      
      
      
      
      it = it+1;
      
      if(verbosity > 5)
        figure(23);
        plot(s0,eta);
        hold on;
        figure(24);
        plot(R,Z,'g');
        hold on;
        pause;
      end
  end
  
 
 
 

 
 if(RES>10)
     disp(' divergence !!!! ' )
     RES
     it
     conv = 0;
     R = R0;
     Z = Z0;
     P = P0
     V = V0
 else
 conv = 1;    
 end
     
    
 %else 
    V = 0;
    for i=2:N
    V = V+ ( (R(i)+R(i-1))^2/4 + (R(i)-R(i-1))^2/12 )*(Z(i)-Z(i-1))*pi;
    end;
%   V
%   P
%   it
%   RES

if(conv==1)
   figure(11);
   plot(R,Z,'r-');
   hold on; 
   conv = 1;
   
   figure(21);
   hold on;
   plot(P,V/Vref,'r.');

figure(31);   
   hold on;
   if(discretization=='FD')
     plot(P,real(lambdaminP),'xr',P,real(lambdaminV),'+b')
   else
     plot(P,real(lambdaminP),'or',P,real(lambdaminV),'ob')    
   end
   
end


 
   
 
 
