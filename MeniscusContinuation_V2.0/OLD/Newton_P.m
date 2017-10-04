function [R,Z,P,V,conv,lambdaminV,lambdaminP,lambdamin1] = Newton_P(R0,Z0,P0,V0,NewP)
% computation of equilibrium solution through Newton method
% input : 
% [R0,Z0] "guess" to start the iteration
% Bo Bond number ( radius divided by capilary length )
% P : nondimensional pressure at nozzle (height of liquid column divided by capilary length)
% output : 
% [R,Z] curve for equilibrium solution
% V  volume

global Bo RR dpdz alphac Weta Wxi verbosity beta discretization Vref; 

P = NewP;

N = length(R0);

%  if(R0(N)~=Bo) % rescaling ofuess so that radius of nozzle = Bo, if not the case 
%    Z0 = Z0*Bo/R0(N);
%    R0 = R0*Bo/R0(N);
%  end

  R = R0;
  Z = Z0;
  RES = 1;
  epsilon = 1e-6;
  it = 0;
  itmax = 100;
  
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
      [eta,RES,lambdamin,lambdamax] = calceta_P(R,Z,alpha,ds,s0,Ka,Kb,dpdz,P);

     
    %methode 1 : eta puis recalcule geometrie puis xi
    if(verbosity < 10)
      R = R + Weta*eta.*sin(alpha);
      Z = Z - Weta*eta.*cos(alpha);
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
      Z = Z  + Wxi*xi.*sin(alpha)-ZN;
     
    else % same with plots for debug
            figure(33);
            plot(R,Z,'r+-')
            hold on;
            R = R + Weta*eta.*sin(alpha);
            Z = Z - Weta*eta.*cos(alpha);
            figure(33);
            plot(R,Z,'bo-')
            if (beta<0)
                xiN = 0;
            else
                [alpha,ds,s0,Ka,Kb] = calcgeom(R,Z);
                xiN = -Z(N)/sin(alpha(N));
            end
            xi = calcxi(R,Z,0,xiN);
            R = R  + Wxi*xi.*cos(alpha);
            Z = Z  + Wxi*xi.*sin(alpha);
            figure(33);
            plot(R,Z,'go-')
    end

       % methode 2 : eta et xi en simultane
       %%%% MARCHE PAS !
%figure(33);
%    plot(R,Z,'r+-')
%    hold on;
    
%      R = R + Weta*eta.*sin(alpha);
%      Z = Z - Weta*eta.*cos(alpha);
%     figure(33);
%      plot(R,Z,'bo-')
    
%       if (beta<0)
%          xiN = 0;
%       else
%          xiN = eta(N)/tan(alpha(N))
%      end
%      xi = calcxi(R,Z,0,xiN)

%      R = R  + Wxi*xi.*cos(alpha);
%      Z = Z  + Wxi*xi.*sin(alpha);
%      figure(33);
%       plot(R,Z,'go-')
%      pause ;


      
      
      if(verbosity > 10)
          RES
      figure(435);
       hold on;
 %      [alpha,ds,s0,Ka,Kb] = calcgeom(R,Z);
       ds(N) = ds(N-1);
       plot(s0,ds-s0(N)/N,'b-x')
      s0apres = s0(N);
       
       figure(436); hold on;
       plot(s0,eta,s0,xi);
        pause;
      end
       
    
      

        if(verbosity > 5)
        figure(23);
        plot(s0,eta);
        hold on;
         
        figure(24);
        plot(R,Z,'g');
        hold on;
        
        figure(25);
        plot(s0,alpha,'b-+');
        hold on;
        
        figure(26);
        plot(s0,Ka,s0,Kb);
        hold on;
        
       
        end
        
         [alpha,ds,s0,Ka,Kb] = calcgeom(R,Z);
         it = it+1
  end
  
 if(RES>10)
     disp(' divergence !!!! ')
     it
     RES
     conv = 0;
     R = R0;
     Z = Z0;
     P = P0;
     V = V0;
     lambdaminP = 0;
     lambdaminV = 0;
     lambdamin1 = 0;
     
 else
     conv = 1;
     V = 0;
     for i=2:N
     V = V+ ( (R(i)+R(i-1))^2/4 + (R(i)-R(i-1))^2/12 )*(Z(i)-Z(i-1))*pi;
     end;
     
     if (verbosity > 1)
     P
     V
     end
 end
 

    
% else 
  
   
   %P
   %V
   %it
   %RES
   %lambdamin
   %lambdamax
   
   if(conv==1)
   figure(11);
   plot(R,Z,'r-');
   hold on;   
   
   figure(20);
   hold on;
   plot(P,V,'r.');
       
   
   [eta,P1,RES,lambdaminV,lambdamaxV,lambdaminP,lambdamaxP,lambdamin1,lambdamax1] = calceta_V(R,Z,P,alpha,ds,s0,Ka,Kb,dpdz,V);
%   eta 
%   P1
   
 figure(31);   
   hold on;
   if(discretization=='FD')
     plot(P,real(lambdaminP),'xr',P,real(lambdamaxP),'or',  ...
          P,real(lambdaminV),'+b',P,real(lambdamaxV),'ob',  ...
          P,real(lambdamin1),'+g',P,real(lambdamax1),'og'  )
   else
     plot(P,real(lambdaminP),'xr',P,real(lambdamaxP),'or',  ...
          P,real(lambdaminV),'+b',P,real(lambdamaxV),'ob',  ...
          P,real(lambdamin1),'+g',P,real(lambdamax1),'og'  )
   end
   
   figure(32);   
   hold on;
   if(discretization=='FD')
     plot(V/Vref,real(lambdaminP),'xr',V/Vref,real(lambdamaxP),'or',  ...
          V/Vref,real(lambdaminV),'+b',V/Vref,real(lambdamaxV),'ob',  ...
          V/Vref,real(lambdamin1),'+g',V/Vref,real(lambdamax1),'og'  )
   else
     plot(V/Vref,real(lambdaminP),'xr',V/Vref,real(lambdamaxP),'or',  ...
          V/Vref,real(lambdaminV),'+b',V/Vref,real(lambdamaxV),'ob',  ...
          V/Vref,real(lambdamin1),'+g',V/Vref,real(lambdamax1),'og'  )
   end
   
   end
    
  