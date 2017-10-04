
classdef meniscus
    
    properties
        % physical properties
        nbdim = 3
        dpdz = 1
        Bo % Bond number ; for drops/bubbles
        L  % Length ; for bridges
        typestart = 'axis' 
        typeend = 'pinned'
        beta = -1
        % number of points 
        N
        % main properties of the meniscus : geometry, pressure, volume
        R
        Z
        P 
        V
        % geometrical properties : angle, curvature, curvilinear abcidia
        alpha 
        Ka
        Kb
        s0
        ds
        % matrices
        AV
        A1
        YoungLap
        % deformations in the iteration processes
        eta
        xi
        xiZ
        % numerical parameters 
        discretization = 'FD'
        paramscan = 'Bo'
        conv = 1;
        verbosity = 0;
        Weta = 1;
        Wxi = 1;
    end
    
    methods
        function m = meniscus(N,Rdisk) % constructor for a plane, circular interface
            m.N = N;
            t = 0:1/N:1;
            m.R = Rdisk*t;
            m.Z = 0*t;
            m.V = 0;
            m.P = 0;
        end % function meniscus
        
         function m = adapt_P(mans,Pnew);
            m = mans;
            m = calcgeom(m);
            m.conv = 0;
            m.P = Pnew;
            N = m.N;
  
  RES = 1;
  epsilon = 1e-6;
  it = 0;
  itmax = 100;
  
   if(verbosity > 2)
    figure(24);plot(R,Z,'r');hold on;
   end
   
  
  while ((RES > epsilon)&&(it<itmax)&&(RES<100))
      % Newton iteration loop
      m = calcgeom(m);
      if(max(abs(m.ds))>1)
          RES = 20000;
      end
      m = calcmat(m);
      m.eta = m.A(1:N,1:N)\m.YoungLap;

%    if(verbosity < 10)
      m.R = m.R + m.Weta*m.eta.*sin(m.alpha);
      m.R = m.Z - m.Weta*m.eta.*cos(m.alpha);
       if (beta<0)
          xi0 = 0; xiN = 0;
          ZN = 0;
       else
          m = calcgeom(m);
%          xiN = -Z(N)/sin(alpha(N)); ZN = 0; % methode 1 : marche pour 45<beta<135
          xiN = 0; xiN = 0;  ZN = Z(N); %methode 2
      end
      m = calcxi(m,xi0,xiN);
      m.R = m.R  + m.Wxi*m.xi.*cos(m.alpha);
      m.Z = m.Z  + m.Wxi*m.xi.*sin(m.alpha)-ZN;
     
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

      if(verbosity > 10)
          RES
          figure(435); hold on; plot(s0,ds-s0(N)/N,'b-x')
          figure(436); hold on;plot(s0,eta,s0,xi);
          pause;
      end
       
    
      

        if(verbosity > 5)
        figure(23);plot(m.s0,m.eta);hold on;
        figure(24);plot(m.R,m.Z,'g');hold on;
        figure(25);plot(m.s0,m.alpha,'b-+');hold on;
        figure(26);plot(m.s0,m.Ka,m.s0,m.Kb);hold on;      
        end
        
         m = calcgeom(m);
         it = it+1;
  end % end of Newton iteration loop
  
 if(RES>10)
     fprintf(' Divergence in Newton P ! \n')
     fprintf(' Previous (P,V) : ( %f , %f ) ; expected new P : %f \n\n',m.P,m.V,Pnew); 
     it
     RES
     m = mans;
     m.conv = 0;
 else
     m.conv = 1;
     if (verbosity > 0)
     m.P
     m.V
     end
 end
 

   
   if(m.conv==1)
   figure(11);plot(m.R,m.Z,'r-');hold on;   
   
   figure(20);hold on;plot(m.P,m.V,'r.');
   end
         end
            
            
           
        
           
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
    
             m.alpha(N) =2*m.alpha(N-1)-m.alpha(N-2)
            %m.alpha(N) = m.alpha(N-2)+2*m.ds(N-1)*(m.Ka(N-1)+m.Kb(N-1)-sin(m.alpha(N-1))/m.R(N-1));
            m.Ka(N) =2*m.Ka(N-1)-m.Ka(N-2); 
            m.Kb(N) =2*m.Kb(N-1)-m.Kb(N-2); 
    
            if(abs(m.R(1))>1e-6)
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
                m.V = V+ ( (R(i)+R(i-1))^2/4 + (R(i)-R(i-1))^2/12 )*(Z(i)-Z(i-1))*pi;
            end;
        
        end % function calcgeom
        
        
        function m = calcmat(m)
        % computation of matrices Ap, Av, A1 and LHS
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
            %if (i~=N-1)
            A(i,i+1) = 2*ds(i-1)/(ds(i)*ds(i-1)*(ds(i)+ds(i-1)));
            %end
            % T0r/r* d eta / d s0
            if (m.nbdim==3)
                A(i,i-1) = A(i,i-1) -ds(i)^2/((ds(i)+ds(i-1))*ds(i)*ds(i-1))*cos(alpha(i))/R(i);
                A(i,i) = A(i,i) +  (ds(i)^2-ds(i-1)^2)/((ds(i)+ds(i-1))*ds(i)*ds(i-1))*cos(alpha(i))/R(i);
                A(i,i+1) = A(i,i+1) +ds(i-1)^2/((ds(i)+ds(i-1))*ds(i)*ds(i-1))*cos(alpha(i))/R(i);
            end
            %  (Ka^2 + Kb^2) eta
            A(i,i) = A(i,i) + 1*(Ka(i)^2+Kb(i)^2);
            % + N0z*eta*dpdz%
            A(i,i) = A(i,i) + dpdz*cos(alpha(i));
        end;
        
        % Right-hand side
        for i=2:N-1
            RHS(i) = Ka(i)+Kb(i)-P+dpdz*Z(i);
        end
        
        
        % case i=1
        if (abs(R(1)) < 1e-6) %for bubble/drop : point 1 is on symmetry axis
            if (nbdim==3)
                A(1,1) = -4/(ds(1)^2);  %%% facteur 4 car 1/R d eta / ds = d^2 eta / d s^2
                A(1,2) = 4/(ds(1)^2);
            else
                A(1,1) = -2/(ds(1)^2);  %%% facteur 4 car 1/R d eta / ds = d^2 eta / d s^2
                A(1,2) = 2/(ds(1)^2);
            end
            A(1,1) = A(1,1)+(Ka(1)^2+Kb(1)^2)+dpdz;
            RHS(1) = Ka(1)+Kb(1)-P+dpdz*Z(1);
            R(1) = 0; % in case of trouble...
        else % for bridge : point 1 is pinned
            A(1,1) = -1;
            A(1,2) = 0;
            RHS(1) = 0;
        end
        
        
        if (beta<0)
            % "pinned" case
            A(N,N) = -1;
            RHS(N) = 0;
        else
            A(N,N-2) =  1/(2*ds(N-1));
            A(N,N-1) = -4/(2*ds(N-1));
            A(N,N)   =  3/(2*ds(N-1));
            
            %A(N,N-1) = -1/(ds(N-1));
            %A(N,N)   =  1/(ds(N-1));
            
            RHS(N)   =  (alpha(N)-beta);
        end
        
        % last equation : volume variation
        % WARNING works only in 3D !
        
        for i=2:N-1
            A(N+1,i) =  (ds(i-1)/2*(2*R(i)+R(i-1))/3 + ds(i)/2*(2*R(i)+R(i+1))/3)*2*pi*TGVvol;
        end;
        %if (R(1)~=0)
        A(N+1,1) = ds(1)/2*(2*R(1)+R(2))/3*2*pi*TGVvol;
        %end
        A(N+1,N) = ds(N-1)/2*(2*R(N)+R(N-1))/3*2*pi*TGVvol;
        
        % last column : effet of dP
        for i=2:N-1
            A(i,N+1) = +1;
        end;
        if(R(1)==0)
            A(1,N+1) = +1;
        end
        
        m.Av = A;
        m.YoungLap = RHS;
        
        
        
        AA1 = A(1:N,1:N);
        for i=2:N
            AA1(i,i) = AA1(i,i) - 1/R(i)^2;
        end
        AA1(1,1) = 1; AA1(1,2) = 0;
        AA1(N,N) = 1; AA1(N,N-1) = 0;
        m.A1 = AA1;
        end
        
       function m = calcxi(m,xi0,xiN);
            m = calcgeom(m);
            N = m.N;
            curvelength = m.s0(N);
            for i=1:N
                m.xi(i)=xi0+(i-1)/(N-1)*(curvelength+xiN)-m.s0(i);
            end
       end       
    end
end