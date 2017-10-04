
classdef meniscus
    
    properties
        %%%% physical properties
        nbdim = 3
        dpdz = 1
        Bo % Bond number ; for drops/bubbles
        L  % Length ; for bridges
        typestart = 'axisX' % possible types : 'axisX', 'pined', 'angle' 
        typeend = 'pined'
        beta = 0
        %%%% number of points 
        N
        %%%% main properties of the meniscus : geometry, pressure, volume
        R
        Z
        P 
        V
        Vans
        Pans
        Vtarget
        %%%% geometrical properties : angle, curvature, curvilinear abcidia
        alpha 
        Ka
        Kb
        s0
        ds
        %%%% matrices
        AV
        AP
        AS
        A1
        YoungLap
        RHS % right-hand-side for newton iteration
        %%% Eigenvalues (for stability criteria)
        lambdaminV
        lambdaminP
        lambdamin1
        %%%% deformations in the iteration processes
        eta
        dP
        xi
        xiZ
        %%%% numerical parameters
        discretization = 'FD';
        paramscan = 'Bo';
        conv = 1;
        verbosity = 5;
        Weta = 1;
        Wxi = 1;
        NewtonEuler = 0;
        typecont = 'P'; 
        % variable for continuation : 'P' (pressure), "V" (volume) or "S" (arclength)
        % Additional parameters for arclength continuation :
        dS;
        anglePV=0;
        dVdS;
        dPdS;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% constructor for a flat meniscus
        function m = meniscus(type,N,param) 
            % type = flat => param is radius
            % type = cylinder => param is length
            
            if (strcmpi(type,'flat')==1)
                m.N = N;
                R0 = param(1);
                t = 0:1/(N-1):1;
                m.R = param*t;
                m.Z = 0*t;
                m.V = 0;
                m.P = 0;
                m.typestart = 'axisX';
                m.typeend = 'pined';
                m.dPdS = 1;
                m.dVdS = 1;
                m = calcgeom(m);
                m = calcmat(m);
                m = calcNewtonEuler(m);
                close all;
                plotmeniscus(m,10,'k'); 
            
            elseif (strcmpi(type,'cylinder')==1)
                L = param(1);
                m.N = N;
                t = 0:1/(N-1):1;
                m.R = ones(size(t));
                m.Z = t*L;
                m.V = pi*L;
                m.P = 1;
                m.typestart = 'pined';
                m.typeend = 'pined';
                m.dpdz = 0;
                m.dPdS = 1;
                m.dVdS = 1;
                m = calcgeom(m);
                m = calcmat(m);
                m = calcNewtonEuler(m);
                close all;
                plotmeniscus(m,10,'k'); 
            
            elseif (strcmpi(type,'sphere')==1)
                R0 = param(1);
                m.beta = param(2);
                m.N = N;
                t = 0:1/(N-1):1;
                t = (m.beta)*t;
                m.R = R0*sin(t);
                m.Z = R0*(cos(t(N))-cos(t));
                m.P = 2/R0;
                m.typestart = 'axisX';
                m.typeend = 'angle';
                m.dpdz = 1*0; %%%%% REMETTRE 1 APRES
                m.dPdS = 1;
                m.dVdS = 1;
                m = calcgeom(m);
                m = calcmat(m);
%                m = calcNewtonEuler(m); % Remettre en route ensuite
                close all;
                plotmeniscus(m,10,'k'); 
            else
                 disp(' Meniscus construction : invalid type !');
            end    
           m.lambdaminV = -1;
           m.lambdaminP = -1;
           m.lambdamin1 = -1;
        end % function meniscus
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function plotmeniscus(m,numfig,plotstyle);        
            figure(numfig);
            plot(m.R,m.Z,plotstyle);
            hold on;
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function m = loop(m,type,dX,Nstep);
              m = step(m,type,0);
              Vtab = [m.V];
              Ptab = [m.P];
              itloop = 0;
              m.typecont=type;
              while (itloop<Nstep && m.conv==1)
                m = step(m,type,dX);
                m = stab(m);
                itloop = itloop+1;
                if(m.conv==1)
                    Vtab = [Vtab m.V];
                    Ptab = [Ptab m.P];
                end   
              end
            
            figure(21);
            hold on;
            plot(Ptab,Vtab,'r');
            plotmeniscus(m,10,'k');
            m.conv=1;
        end
        


         function m = step(mans,type,dX);
            m = mans;
            m = calcgeom(m);
            m.conv = 0;
            m.typecont = type;
            m.Vans = m.V;
            m.Pans = m.P;
%            Pans = m.P;
%            m.P = Pnew;
            N = m.N;
           
  
  RES = 1;
  epsilon = 1e-6;
  it = 0;
  itmax = 100;
  
   if(m.verbosity > 2)
      plotmeniscus(m,24,'r');
   end
   
%  if (m.NewtonEuler==1) %% a revoir
%      m = calcNewtonEuler(m);
%        m.R = m.R + dX*m.Weta*m.eta.*sin(m.alpha);
%        m.Z = m.Z - dX*m.Weta*m.eta.*cos(m.alpha);
%        m.P = m.P+dX*dP;
%  end
  
   if (m.typecont == 'P')
   m.P = m.P+dX; 
   end
   if (m.typecont == 'V')
   m.Vtarget = m.V+dX; 
   end
   if (m.typecont == 'S')
        m.dS = dX; 
        m = calcmat(m);
        m = calcNewtonEuler(m);
        m.R = m.R + m.dS*m.dPdS*m.eta.*sin(m.alpha);
        m.Z = m.Z - m.dS*m.dPdS*m.eta.*cos(m.alpha);
        m.P = m.P+m.dS*m.dPdS;
        Vans = m.V;
        m = calcgeom(m);
        %debug
        %figure(20);
        %plot(m.P,m.V,'bx');
   end
   
  while ((RES > epsilon)&&(it<itmax)&&(RES<100))
      % Newton iteration loop
      
      m = calcgeom(m);
      if(max(abs(m.ds))>1)
          RES = 20000;
      end
      m = calcmat(m);
      plotmeniscus(m,55,'k+-'); % a virer
      if (m.typecont == 'P')
      m = calcRHS(m);
      % ancienne methode
      m.eta = (m.AP(1:m.N,1:m.N)\m.RHS(1:m.N)')';
      
      figure(55);
      plot(m.R + m.Weta*m.eta.*sin(m.alpha),m.Z - m.Weta*m.eta.*cos(m.alpha),'g+-');
      pause;
      figure(436); hold on;plot(m.s0,m.eta,'g+');
      
   %      m.P = m.P+dP;

% Nouvelle methode    
      
      etaXX = (m.AP(1:2*m.N,1:2*m.N)\m.RHS(1:2*m.N)')';
      m.eta = etaXX(1:m.N);
      m.xi = etaXX(m.N+1 : 2*m.N); 
      dP = 0;
      RES = max(abs(m.RHS(1:m.N)));
      end
      
      if (m.typecont == 'V')
      m = calcRHS(m);
      etaXX = (m.AV(1:m.N+1,1:m.N+1)\m.RHS(1:m.N+1)')';
      m.eta = etaXX(1:m.N);
      dP = etaXX(m.N+1);
      RES = max(abs(m.RHS(1:m.N)));
      end
      
      if (m.typecont == 'S')
      m = calcRHS(m);
      etaXX = (m.AS(1:m.N+1,1:m.N+1)\m.RHS(1:m.N+1)')';
      m.eta = etaXX(1:m.N);
      dP = etaXX(m.N+1);
      RES = max(abs(m.RHS(1:m.N)));
      end
      
%      figure(333);
%      plot(m.s0(1:m.N),m.RHS(1:m.N),'b+',m.s0(1:m.N),m.eta(1:m.N),'r');
%      hold on;


     figure(436); hold on;plot(m.s0,m.eta,'r+');
     pause;
     plot(m.s0,m.xi,'b+');
     
      m.R = m.R + m.Weta*m.eta.*sin(m.alpha);
      m.Z = m.Z - m.Weta*m.eta.*cos(m.alpha);
      m.P = m.P+dP;

% correction to ensure equal spacing of gridpoints      
      
      plotmeniscus(m,55,'r+');
      pause ;
      
 %     m = calcxi(m);
 %     ZN = m.Z(m.N);
 %     m.R = m.R  + m.Wxi*m.xi.*cos(m.alpha);
 %     m.Z = m.Z  + m.Wxi*m.xi.*sin(m.alpha)-ZN;
 
 %           xi0 = 0; xiN = 0;  
 %           m = calcgeom(m);
            N = m.N;
 %           curvelength = m.s0(m.N);
 %           for i=1:N
 %               m.xi(i)=xi0+(i-1)/(N-1)*(curvelength+xiN)-m.s0(i);
 %           end
            
       m.R = m.R  + m.Wxi*m.xi.*cos(m.alpha);
       m.Z = m.Z  + m.Wxi*m.xi.*sin(m.alpha);     
       plotmeniscus(m,55,'b+');
       pause ;
       
       
      
      %debug
%      figure(20);
%      plot(m.P,m.V,'g+');
%    else % same with plots for debug
%            figure(33);
%            plot(R,Z,'r+-')
%            hold on;
%            R = R + Weta*eta.*sin(alpha);
%            Z = Z - Weta*eta.*cos(alpha);
%            figure(33);
%            plot(R,Z,'bo-')
%            if (beta<0)
%                xiN = 0;
%            else
%                [alpha,ds,s0,Ka,Kb] = calcgeom(R,Z);
%                xiN = -Z(N)/sin(alpha(N));
%            end
%            xi = calcxi(R,Z,0,xiN);
%            R = R  + Wxi*xi.*cos(alpha);
%            Z = Z  + Wxi*xi.*sin(alpha);
%            figure(33);
%            plot(R,Z,'go-')
%    end

      if(m.verbosity > 10)
          RES
          figure(435); hold on; plot(m.s0,m.ds-m.s0(N)/N,'b-x')
          figure(436); hold on;plot(m.s0,m.eta,m.s0,m.xi);
          pause;
      end
       
    
      

        if(m.verbosity > 5)
        figure(23);plot(m.s0,m.eta);hold on;
        figure(24);plot(m.R,m.Z,'g');hold on;
        figure(25);plot(m.s0,m.alpha,'b-+');hold on;
        figure(26);plot(m.s0,m.Ka,m.s0,m.Kb);hold on;      
        end
        
         m = calcgeom(m);
         it = it+1;
  end % end of Newton iteration loop
       
  
  
%  if (m.typecont == 'P')
%  m.P = Pnew;  
%  end
%  if (m.typecont == 'V')
%  m.V = Vnew;  
%  end
  
  
 if(RES>10)
     fprintf(' Divergence in Newton iteration ! \n')
     m = mans;
     m.conv = 0;
     fprintf(' Previous (P,V) : ( %f , %f ) \n\n',m.P,m.V); 
 else
     m.conv = 1;
     if (m.verbosity > 0)
         fprintf(' New converged meniscus shape for (P,V) = ( %f , %f ); it = %i \n',m.P,m.V,it); 
     end
 end
 

   
   if(m.conv==1)
   figure(11);plot(m.R,m.Z,'r-');hold on;   
   
   figure(20);hold on;plot(m.P,m.V,'r.');
   end
         end
  
         
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
        
           
        function m = calcgeom(m);
            m.ds(1) = sqrt((m.R(2)-m.R(1))^2+(m.Z(2)-m.Z(1))^2);
            m.s0(1) = 0;
            N = m.N;
            
            % curvilinear abcidia
            for i=2:N-1
                m.ds(i) = sqrt( (m.R(i+1)-m.R(i))^2+(m.Z(i+1)-m.Z(i))^2);
                m.s0(i) = m.s0(i-1)+m.ds(i-1);
            end
            m.s0(N) = m.s0(N-1)+m.ds(N-1);
            m.ds(N) = m.ds(N-1); % not significant but to avoid possible bugs
            
            % angle at center of sides (intermediate array)
            for i=1:N-1
                alphaS(i) =  atan2(m.Z(i+1)-m.Z(i),m.R(i+1)-m.R(i)); % angle at the middle of the side
            end       
            unwrap(alphaS); 
            %   angle and curvature at nodes
            for i=2:N-1
                m.alpha(i) =  (m.ds(i)*alphaS(i-1) +  m.ds(i-1)*alphaS(i))/(m.ds(i-1)+m.ds(i));
                m.Ka(i) = (alphaS(i)-alphaS(i-1))/((m.ds(i)+m.ds(i-1))/2);
            if (m.nbdim == 3) 
                m.Kb(i) = sin(m.alpha(i))/m.R(i);
            else
                m.Kb(i) = 0;
            end
            end;
    
            % i = N : pinned point (geom. quantities only needed for visualization and for angle-fixed case)
    
             m.alpha(N) =2*m.alpha(N-1)-m.alpha(N-2);
            %m.alpha(N) = m.alpha(N-2)+2*m.ds(N-1)*(m.Ka(N-1)+m.Kb(N-1)-sin(m.alpha(N-1))/m.R(N-1));
            m.Ka(N) =2*m.Ka(N-1)-m.Ka(N-2); 
            m.Kb(N) =2*m.Kb(N-1)-m.Kb(N-2); 
    
            if(m.typestart=='pined')
                %  i = 1 : pinned point
                m.alpha(1) = 2*m.alpha(2)-m.alpha(3);
                m.Ka(1) = 2*m.Ka(2)-m.Ka(3);
                m.Kb(1) = 2*m.Kb(2)-m.Kb(3);
            else
                %  i = 1 : axis
                m.alpha(1) = 0.;
                m.Ka(1) = 2*atan2(m.Z(2)-m.Z(1),m.R(2)-m.R(1) )/m.ds(1); % facteur 2 rajouté
                if (m.nbdim == 3)
                    m.Kb(1) = m.Ka(1);
                end
            end
            unwrap(m.alpha);
            
            m.V = 0;
            for i=2:N
                m.V = m.V+ ( (m.R(i)+m.R(i-1))^2/4 + (m.R(i)-m.R(i-1))^2/12 )*(m.Z(i)-m.Z(i-1))*pi;
            end;
        
        end % function calcgeom
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function m = calcmat(m)
        % computation of matrices Ap, Av, A1 
        % Equation : -K1(eta)+N0z*eta*dpdz = (K0 - (P+dpdz*z)) 
        % Resolution as a matricial equation : A * eta = RHS

        A = 0;
        N = m.N;
        ds = m.ds;
        R = m.R;
        Z = m.Z;
        alpha = m.alpha;
        Ka = m.Ka;
        Kb = m.Kb;
        
        for i=2:N-1
            % d^2 eta / d^2 s0
            A(i,i-1) = 2*ds(i)/(ds(i)*ds(i-1)*(ds(i)+ds(i-1)));
            A(i,i) = -2/(ds(i)*ds(i-1));
            A(i,i+1) = 2*ds(i-1)/(ds(i)*ds(i-1)*(ds(i)+ds(i-1)));
            % T0r/r* d eta / d s0
            if (m.nbdim==3)
                A(i,i-1) = A(i,i-1) -ds(i)^2/((ds(i)+ds(i-1))*ds(i)*ds(i-1))*cos(alpha(i))/R(i);
                A(i,i) = A(i,i) +  (ds(i)^2-ds(i-1)^2)/((ds(i)+ds(i-1))*ds(i)*ds(i-1))*cos(alpha(i))/R(i);
                A(i,i+1) = A(i,i+1) +ds(i-1)^2/((ds(i)+ds(i-1))*ds(i)*ds(i-1))*cos(alpha(i))/R(i);
            end
            %  (Ka^2 + Kb^2) eta
            A(i,i) = A(i,i) + 1*(Ka(i)^2+Kb(i)^2);
            % + N0z*eta*dpdz%
            A(i,i) = A(i,i) + m.dpdz*cos(alpha(i));
        end;
        

        if (strcmpi(m.typestart,'axisX')==1) %for bubble/drop : point 1 is on symmetry axis
            if (m.nbdim==3)
                A(1,1) = -4/(ds(1)^2);  %%% facteur 4 car 1/R d eta / ds = d^2 eta / d s^2
                A(1,2) = 4/(ds(1)^2);
            else
                A(1,1) = -2/(ds(1)^2);  
                A(1,2) = 2/(ds(1)^2);
            end
            A(1,1) = A(1,1)+(Ka(1)^2+Kb(1)^2)+m.dpdz;
        elseif (strcmpi(m.typestart,'pined')==1) % for bridge : point 1 is pinned
            A(1,1) = -1;
            A(1,2) = 0;
        end
        
        if (strcmpi(m.typeend,'pined')==1)
            % "pined" case
            A(N,N) = -1;
            RHS(N) = 0;
        elseif (strcmpi(m.typeend,'angle')==1)
            A(N,N-2) =  1/(2*ds(N-1));
            A(N,N-1) = -4/(2*ds(N-1));
            A(N,N)   =  3/(2*ds(N-1));
        end
        
         %%% Part of matrix AP corresponding to xi
         % - [ d ( Ka + Kb )/ds0 + dpdz*sin(alpha) ] xi
         for i = 2:N-1
             A(i,N+i) = -(((m.Ka(i+1)+m.Kb(i+1))-(m.Ka(i-1)+m.Kb(i-1)))/(ds(i)+ds(i-1)));
             A(i,N+i) = A(i,N+i) - m.dpdz*sin(m.alpha(i));    % sign ?
         end
        
        %  
        % d^2 xi / d^2 s0
          for i = 2:N-1
            A(N+i,N+i-1) =  2*ds(i)/(ds(i)*ds(i-1)*(ds(i)+ds(i-1)));
            A(N+i,N+i)   = -2/(ds(i)*ds(i-1));
            A(N+i,N+i+1) =  2*ds(i-1)/(ds(i)*ds(i-1)*(ds(i)+ds(i-1)));
         end
         % ( d Ka / d s0 ) eta + Ka (d eta / d s0)
          for i = 2:N-1
                A(N+i,i-1) = Ka(i)*(-ds(i)^2/((ds(i)+ds(i-1))*ds(i)*ds(i-1)));
                A(N+i,i)   = Ka(i)*(ds(i)^2-ds(i-1)^2)/((ds(i)+ds(i-1))*ds(i)*ds(i-1));
                A(N+i,i+1) = Ka(i)*ds(i-1)^2/((ds(i)+ds(i-1))*ds(i)*ds(i-1));
                A(N+i,i)   = A(N+i,i) + (m.Ka(i+1)-m.Ka(i-1))/(ds(i)+ds(i-1));
          end
        
          if (strcmpi(m.typestart,'pined')==1)
            A(N+1,N+1) = 1;
          elseif (strcmpi(m.typestart,'axisX')==1)
            A(N+1,N+1) = 1;
          end
          
          if (strcmpi(m.typeend,'pined')==1)
            A(2*N,2*N) = 1;
          elseif (strcmpi(m.typeend,'angle')==1) 
            A(2*N,N) = 1;
            A(2*N,2*N) = -tan(m.beta);
          end
         
          m.AP = A;
          
        %%% END OF MATRIX AP ; NEXT IS FOR AV or AS
         m.AV(1:N,1:N) = A(1:N,1:N);
         m.AS(1:N,1:N) = A(1:N,1:N);        
       
        % last column  for AV  / AS : effet of dP
        for i=2:N-1
            m.AV(i,N+1) = +1;
        end;
        if(R(1)==0)
            m.AV(1,N+1) = +1;
        end
        m.AS(1:N,N+1) = m.AV(1:N,N+1);
        m.AV(N+1,N+1) = 0;
         
        % last line for AV : Volume increment
        for i=2:N-1
            m.AV(N+1,i) =  (ds(i-1)/2*(2*R(i)+R(i-1))/3 + ds(i)/2*(2*R(i)+R(i+1))/3)*2*pi;
        end;
        m.AV(N+1,1) = ds(1)/2*(2*R(1)+R(2))/3*2*pi;
        m.AV(N+1,N) = ds(N-1)/2*(2*R(N)+R(N-1))/3*2*pi;

        % last line for AS : increment in arclength
%        m = calcNewtonEuler(m);  % mandatory to get dPdS and dVdS
%        m.AS(N+1,1:N) = m.dVdS*m.AV(N+1,1:N);
%        m.AS(N+1,N+1) = m.dPdS;
        
 % Matrix A1 : for non-axisymetric instabilities       
       
        AA1 = A(1:N,1:N);
        for i=2:N
            AA1(i,i) = AA1(i,i) - 1/R(i)^2;
        end
        AA1(1,1) = -1; AA1(1,2) = 0;
        AA1(N,N) = -1; AA1(N,N-1) = 0;
        m.A1 = AA1;
        end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        
    function m = calcRHS(m);
        
        % Laplace/Young operator for central points
        for i=2:m.N-1
            m.RHS(i) = m.Ka(i)+m.Kb(i)-m.P+m.dpdz*m.Z(i);
        end
        % starting point
        if (strcmpi(m.typestart,'axisX')==1)
            % for bubble/drop : point 1 is on symmetry axis
            m.RHS(1) = m.Ka(1)+m.Kb(1)-m.P+m.dpdz*m.Z(1);
         end
         if (strcmpi(m.typestart,'pined')==1)
            % for bridge : point 1 is pinned
            m.RHS(1) = 0;
         end 
         % ending point
         if (strcmpi(m.typeend,'pined')==1)
            m.RHS(m.N) = 0;
         end
         if (strcmpi(m.typeend,'angle')==1)
            m.RHS(m.N)   =  (m.alpha(m.N)-m.beta);
         end
         for i=1:m.N
            m.RHS(m.N+i) = 0;
        end
        
         
         %%%%% A REVOIR ENSUITE 
         % For V continuation : volume constraint
         if (m.typecont=='V')
            m.RHS(m.N+1) = m.Vtarget-m.V;
         end
         
         % For S continuation : arclength step
         if (m.typecont=='S')
%            dSi2 = (m.V-m.Vans)^2+(m.P-m.Pans)^2; 
%            m.RHS(m.N+1) = m.dS^2-dSi2;
            m.RHS(m.N+1) = 0;
         end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function m = calcNewtonEuler(m);       
     % improved guess using NewtonEuler method
     if (m.typecont == 'P'|m.typecont == 'S')
        RHS  = -m.AV(1:m.N,m.N+1);
        m.eta = (m.AP\RHS)';
        dP = 1;
        dV = m.AV(m.N+1,1:m.N)*m.eta';
        dS = sqrt(dP^2+dV^2);
        if (m.typecont == 'S')
             dPdSi = dP/dS;
             dVdSi = dV/dS;
           if(dPdSi*m.dPdS+dVdSi*m.dVdS>0)
              m.dPdS = dPdSi;
              m.dVdS = dVdSi;
           else
              m.dPdS = -dPdSi;
              m.dVdS = -dVdSi;
           end
        end
     m.AS(m.N+1,1:m.N) = m.dVdS*m.AV(m.N+1,1:m.N);
        m.AS(m.N+1,m.N+1) = m.dPdS;   
     end
     if (m.typecont == 'V')
        RHS  = [ones(m.N,1)*0;1];
        etaXX =  (m.AV\RHS)';
        m.eta = etaXX(1:m.N);
        dP = etaXX(m.N+1);
        dV = 1;
        dS = sqrt(dP^2+dV^2);
        m.dPdS = dP/dS;
        m.dVdS = dV/dS;
     end
    end
    
     
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       function m = calcxi(m);
          
       if (m.typeend=='pined')
          xi0 = 0; xiN = 0;
       else
          m = calcgeom(m);
          xi0 = 0;
%          xiN = -m.Z(m.N)/sin(m.alpha(m.N)); ZN = 0; % methode 1 : marche pour 45<beta<135
          xi0 = 0; xiN = 0;   %methode 2
      end
            m = calcgeom(m);
            N = m.N;
            curvelength = m.s0(m.N);
            for i=1:N
                m.xi(i)=xi0+(i-1)/(N-1)*(curvelength+xiN)-m.s0(i);
            end
            
       end
       
       
    function m = loopANS(m,typeloop,dX,Nstep)
        m.conv = 1;
        mans = m;
%        Pstart = m.P;
%        PP = Pstart;

        Vtab = [m.V];
        Ptab = [m.P];
        itloop = 0;
        m.typecont=typeloop;
        while ((itloop<=Nstep)&&(m.conv==1))
    
%lambdaminVans = lambdaminV;
%lambdaminPans = lambdaminP;
%lambdamin1ans = lambdamin1;
%Vans = V;
%Pans = P;
    
%[R,Z,P,V,conv,lambdaminV,lambdaminP,lambdamin1] = Newton_P(R,Z,P,V,PP);
        m = step(m,dX);


        if(m.conv==1)
            Vtab = [Vtab m.V];
            Ptab = [Ptab m.P];

%if (lambdaminV*lambdaminVans<0&&abs(lambdaminV)~=1&&abs(lambdaminVans)~=1) % detection bifurcation V
%    VbV = (V*lambdaminVans - Vans*lambdaminV)/(lambdaminVans-lambdaminV);
%    PbV = (P*lambdaminVans - Pans*lambdaminV)/(lambdaminVans-lambdaminV);
%    figure(21);
%    hold on;
%    plot(PbV,VbV/Vref,'+r','MarkerSize',10);
%    figure(10);
%    hold on;
%    plot(R,Z*sign(dpdz),'b')
%    disp('Bifurcation V detectee !')
%    PbV
%    VbV
%    tabV = [tabV,[Bo ; PbV ; VbV]];
%    conv = 0 % (pour diagrammes V/Bo ; a virer ensuite)
%end

%if (lambdaminP*lambdaminPans<0&&abs(lambdaminP)~=1&&abs(lambdaminPans)~=1) % detection bifurcation V
%    VbP = (V*lambdaminPans - Vans*lambdaminP)/(lambdaminPans-lambdaminP);
%    PbP = (P*lambdaminPans - Pans*lambdaminP)/(lambdaminPans-lambdaminP);
%    figure(21);
%    hold on;
%    plot(PbP,VbP/Vref,'xr','MarkerSize',10); 
%    figure(10);
%    hold on;
%    plot(R,Z*sign(dpdz),'k:')
%    disp('Bifurcation P detectee !')
%    PbP
%    VbP
%     tabP = [tabP,[Bo ; PbP ; VbP]];
%end

%if (lambdamin1*lambdamin1ans<0&&abs(lambdamin1)~=1&&abs(lambdamin1ans)~=1) % detection bifurcation V
%    Vb1 = (V*lambdamin1ans - Vans*lambdamin1)/(lambdamin1ans-lambdamin1);
%    Pb1 = (P*lambdamin1ans - Pans*lambdamin1)/(lambdamin1ans-lambdamin1);
%    figure(21);
%    hold on;
%    plot(Pb1,Vb1/Vref,'*k','MarkerSize',10); 
%    figure(10);
%    hold on;
%    plot(R,Z*sign(dpdz),'m')
%    disp('Bifurcation m=1 detectee !')
%    Pb1
%    Vb1
%    tab1 = [tab1,[Bo ; Pb1 ; Vb1]];
%end

end

%PP = PP+dP;
itloop = itloop+1;
end



figure(21);
hold on;
plot(Ptab,Vtab,'r');

figure(10);
hold on;
plot(m.R,m.Z,'r');
    end
       
       
       
       
   

    function m = stab(m);
    
          
    lambdaminVans = m.lambdaminV;
    lambdaminPans = m.lambdaminP;
    lambdamin1ans = m.lambdamin1;
        
     EE = eig(m.AP);
     EV = EE(abs(EE)==min(abs(EE)));
     m.lambdaminP = EV(1);
     
     EE = eig(m.AV);
     EV = EE(abs(EE)==min(abs(EE)));
     m.lambdaminV = EV(1);
     
     EE = eig(m.A1);
     EV = EE(abs(EE)==min(abs(EE)));
     m.lambdamin1 = EV(1); 
     
     if(m.verbosity>1) 
        figure(22);
        hold on;
        plot(m.P,real(m.lambdaminP),'r+',m.P,real(m.lambdaminV),'b+',m.P,real(m.lambdamin1),'g+');
        plot(m.P,imag(m.lambdaminP),'rx',m.P,imag(m.lambdaminV),'bx',m.P,imag(m.lambdamin1),'gx');
     end
 
     if (m.lambdaminV*lambdaminVans<0&&abs(m.lambdaminV-lambdaminVans)<1) % detection bifurcation V
     VbV = (m.V*lambdaminVans - m.Vans*m.lambdaminV)/(lambdaminVans-m.lambdaminV);
     PbV = (m.P*lambdaminVans - m.Pans*m.lambdaminV)/(lambdaminVans-m.lambdaminV);
      if (imag(VbV)==0)
        figure(21);
        hold on;
        plot(PbV,VbV,'+b','MarkerSize',10);
        mbif = step(m,'P',PbV-m.P);
        plotmeniscus(mbif,10,'b');
        fprintf('***\n*** V-bifurcation detected for (P,V) = ( %f , %f ) \n***\n', PbV,VbV); 
      end
    end
    
    if (m.lambdaminP*lambdaminPans<0&&abs(m.lambdaminP-lambdaminPans)<1) % detection bifurcation P
      VbV = (m.V*lambdaminPans - m.Vans*m.lambdaminP)/(lambdaminPans-m.lambdaminP);
      PbV = (m.P*lambdaminPans - m.Pans*m.lambdaminP)/(lambdaminPans-m.lambdaminP);
      if (imag(VbV)==0)
        figure(21);
        hold on;
        plot(PbV,VbV,'xb','MarkerSize',10);
        mbif = step(m,'V',VbV-m.V);
        plotmeniscus(mbif,10,'r');
        fprintf('***\n*** P-bifurcation detected for (P,V) = ( %f , %f ) \n***\n', PbV,VbV);    
      end
    end
    if (m.lambdamin1*lambdamin1ans<0&&abs(m.lambdamin1-lambdamin1ans)<1) % detection bifurcation V
     VbV = (m.V*lambdamin1ans - m.Vans*m.lambdamin1)/(lambdamin1ans-m.lambdamin1);
     PbV = (m.P*lambdamin1ans - m.Pans*m.lambdamin1)/(lambdamin1ans-m.lambdamin1);
     figure(21);
     hold on;
     plot(PbV,VbV,'*b','MarkerSize',10);
     plotmeniscus(m,10,'b');
     if (imag(VbV)==0)
       mbif = step(m,'V',VbV-m.V);
       plotmeniscus(mbif,10,'g');    
     end
     fprintf('***\n*** V-bifurcation detected for (P,V) = ( %f , %f ) \n***\n', PbV,VbV);   
    end

%if (lambdaminP*lambdaminPans<0&&abs(lambdaminP)~=1&&abs(lambdaminPans)~=1) % detection bifurcation V
%    VbP = (V*lambdaminPans - Vans*lambdaminP)/(lambdaminPans-lambdaminP);
%    PbP = (P*lambdaminPans - Pans*lambdaminP)/(lambdaminPans-lambdaminP);
%    figure(21);
%    hold on;
%    plot(PbP,VbP/Vref,'xr','MarkerSize',10); 
%    figure(10);
%    hold on;
%    plot(R,Z*sign(dpdz),'k:')
%    disp('Bifurcation P detectee !')
%    PbP
%    VbP
%     tabP = [tabP,[Bo ; PbP ; VbP]];
%end

%if (lambdamin1*lambdamin1ans<0&&abs(lambdamin1)~=1&&abs(lambdamin1ans)~=1) % detection bifurcation V
%    Vb1 = (V*lambdamin1ans - Vans*lambdamin1)/(lambdamin1ans-lambdamin1);
%    Pb1 = (P*lambdamin1ans - Pans*lambdamin1)/(lambdamin1ans-lambdamin1);
%    figure(21);
%    hold on;
%    plot(Pb1,Vb1/Vref,'*k','MarkerSize',10); 
%    figure(10);
%    hold on;
%    plot(R,Z*sign(dpdz),'m')
%    disp('Bifurcation m=1 detectee !')
%    Pb1
%    Vb1
%    tab1 = [tab1,[Bo ; Pb1 ; Vb1]];
%end

    end


    

    end
    end