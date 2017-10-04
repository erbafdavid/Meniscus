
classdef meniscus
    
    properties
        %%%% physical properties
        nbdim = 3
        dpdz = 1
        Bo % Bond number ; for drops/bubbles
        L  % Length ; for bridges
        typestart = 'axisX' % possible types : 'axisX', 'pined' 
        typeend = 'pined'   % possible types : 'pined', 'angle'
        beta = 0            % contact angle, for angle-controled case
        %%%% number of points
        N
        %%%% main properties of the meniscus : geometry, pressure, volume
        R
        Z
        P
        V
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
        RHS % right-hand-side for newton iteration
        %%% Eigenvalues (for stability criteria)
        lambdaminV
        lambdaminP
        lambdamin1
        %%%% deformations in the iteration processes
        eta
        dP
        xi
        %%%% numerical parameters
        Vans
        Pans
        Vtarget
        discretization = 'FD';
        %paramscan = 'Bo';
        conv = 1;
        verbosity = 5;
        Weta = 1;
        Wxi = 1;
        NewtonEuler = 0;
        typecont = 'P';
        % variable for continuation : 'P' (pressure), "V" (volume) or "S" (arclength)
        % Additional parameters for arclength continuation :
        dS;
        dVdS;
        dPdS;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function m = meniscus(type,N,param)
            % Constructor for a guess Meniscus
            % type = flat => param is Radius
            % type = cylinder => param is [Radius,Length,Bond]
            % type = sphere => spherical cap, angle-controled, param is [Rcap,beta]
            
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
                m.dpdz = 1;
            else
                disp(' Meniscus construction : invalid type !');
            end
            
            m = calcgeom(m);
            m = calcmat(m);
            m.dPdS = 1;
            m.dVdS = 0;
            m = calcNewtonEuler(m);
            close all;
            plotmeniscus(m,10,'k:');
            
            m.lambdaminV = -1;
            m.lambdaminP = -1;
            m.lambdamin1 = -1;
            m.conv = 0;
        end % function meniscus
        
       
        
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        function m = step(mans,type,dX);
            % Newton step from a previous meniscus.
            % continuation type may be 'P', 'V', 'S'
            
            m = mans;
            m = calcgeom(m);
            m.conv = 0;
            m.typecont = type;
            m.Vans = m.V;
            m.Pans = m.P;
            N = m.N;
            
            
            RES = 1;
            epsilon = 1e-6;
            it = 0;
            itmax = 100;
  
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
            end
            
            while ((RES > epsilon)&&(it<itmax)&&(RES<100))
                % Newton iteration loop
                
                m = calcgeom(m);
                if(max(abs(m.ds))>1)
                    RES = 20000;
                end
                m = calcmat(m);
                
                if (m.typecont == 'P')
                    m = calcRHS(m);
                    m.eta = (m.AP(1:m.N,1:m.N)\m.RHS(1:m.N)')';
                    dP = 0;
                    RES = max(abs(m.RHS(1:m.N)));
                end
                
                if (m.typecont == 'V')
                    m = calcRHS(m);
                    etaXX = (m.AV(1:m.N+1,1:m.N+1)\m.RHS(1:m.N+1)')';
                    m.eta = etaXX(1:m.N);
                    dP = etaXX(m.N+1);
                    RES = max(abs(m.RHS(1:m.N+1)));
                end
                
                if (m.typecont == 'S')
                    m = calcRHS(m);
                    etaXX = (m.AS(1:m.N+1,1:m.N+1)\m.RHS(1:m.N+1)')';
                    m.eta = etaXX(1:m.N);
                    dP = etaXX(m.N+1);
                    RES = max(abs(m.RHS(1:m.N)));
                end
                
                m.R = m.R + m.Weta*m.eta.*sin(m.alpha);
                m.Z = m.Z - m.Weta*m.eta.*cos(m.alpha);
                m.P = m.P+dP;
                
                % correction to ensure equal spacing of gridpoints and
                % reset Z(N) = 0.
                if (m.typeend=='angle')
                    ZN = m.Z(N);
                else
                    ZN = 0;
                end
                
                m = calcxi(m);
                m.R = m.R  + m.Wxi*m.xi.*cos(m.alpha);
                m.Z = m.Z  + m.Wxi*m.xi.*sin(m.alpha)-ZN;
            
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
            end
              
            if(m.conv==1)
                plotmeniscus(m,11,'r-');                
                figure(21);hold on;plot(m.P,m.V,'r.');
            end
        end % function step
        
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
         function m = loop(m,type,dX,Nstep);
              % Loop over Newton steps, to draw P/V diagrams. 
              % Including eigenvalue co
              % continuation type may be 'P', 'V', 'S'
            m = step(m,type,0);
            Vtab = [m.V];
            Ptab = [m.P];
            plotmeniscus(m,10,'k');
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
            
            figure(20);
            hold on;
            plot(Ptab,Vtab,'r');
            plotmeniscus(m,10,'k');
            m.conv=1;
            
            % for case where next loop will be S-continuation : go forward
            if (m.typecont == 'P')
                m.dPdS=dX;
                m.dVdS=0; 
            end
            if (m.typecont == 'V')
                m.dVdS=dX;
                m.dPdS=0; 
            end
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function m = calcgeom(m);
            % to compute geometrical properties of meniscus :
            % alpha, Ka, Kb, s0, ds
            
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
            % computation of matrices AP, AV, A1
            %
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
                A(i,i) = A(i,i) + m.dpdz*cos(alpha(i));
            end;
            
            
            % case i=1
            if (strcmpi(m.typestart,'axisX')==1) %for bubble/drop : point 1 is on symmetry axis
                if (m.nbdim==3)
                    A(1,1) = -4/(ds(1)^2);  %%% facteur 4 car 1/R d eta / ds = d^2 eta / d s^2
                    A(1,2) = 4/(ds(1)^2);
                else
                    A(1,1) = -2/(ds(1)^2);
                    A(1,2) = 2/(ds(1)^2);
                end
                A(1,1) = A(1,1)+(Ka(1)^2+Kb(1)^2)+m.dpdz;
            else % for bridge : point 1 is pinned
                A(1,1) = -1;
                A(1,2) = 0;
            end
            
            
            if (strcmpi(m.typeend,'pined')==1)
                % "pinned" case
                A(N,N) = -1;
            elseif(strcmpi(m.typeend,'angle')==1)
                A(N,N-2) =  1/(2*ds(N-1));
                A(N,N-1) = -4/(2*ds(N-1));
                A(N,N)   =  3/(2*ds(N-1));
                for i=1:N-1
                    A(i,N) = A(i,N) - m.dpdz*cos(alpha(N));
                end
            end
            
            m.AP = A;
            m.AV(1:N,1:N) = A;
            m.AS(1:N,1:N) = A;
            
            
            %%% END OF MATRIX AP ; NEXT IS FOR AV or AS
            
            
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
            
            if(strcmpi(m.typeend,'angle')==1)
                m.AV(N+1,N) = m.AV(N+1,N)-pi*R(N)^2*cos(m.alpha(N));
            end
            
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
        
        end % function calcmat
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function m = calcRHS(m);
            % computes the right-hand-side in P,V and S continuations
            
            % Laplace/Young operator for central points
            for i=2:m.N-1
                m.RHS(i) = m.Ka(i)+m.Kb(i)-m.P+m.dpdz*(m.Z(i)-m.Z(m.N));
            end
            % starting point
            if (strcmpi(m.typestart,'axisX')==1)
                % for bubble/drop : point 1 is on symmetry axis
                m.RHS(1) = m.Ka(1)+m.Kb(1)-m.P+m.dpdz*(m.Z(1)-m.Z(m.N));
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
            % Guess for S-iteration using NewtonEuler method
            
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
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function m = calcxi(m);
            m = calcgeom(m);
            N = m.N;
            curvelength = m.s0(N);
            for i=1:N
                m.xi(i)=(i-1)/(N-1)*(curvelength)-m.s0(i);
            end
        end
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
             
        function m = stab(m);
            % computes the minimun eigenvalues of AP,AV,A1 and checks if a
            % bifurcation is encountered in the loop process
            
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
            
            if (m.lambdaminV*lambdaminVans<0&&abs(m.lambdaminV-lambdaminVans)<1&&abs(imag(m.lambdaminV))<0.1) % detection bifurcation V
                VbV = (m.V*lambdaminVans - m.Vans*m.lambdaminV)/(lambdaminVans-m.lambdaminV);
                PbV = (m.P*lambdaminVans - m.Pans*m.lambdaminV)/(lambdaminVans-m.lambdaminV);
                    figure(21);
                    hold on;
                    plot(PbV,VbV,'+b','MarkerSize',10);
                if (imag(VbV)==0)
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
                fprintf('***\n*** NonAxi-bifurcation detected for (P,V) = ( %f , %f ) \n***\n', PbV,VbV);
            end
            
          
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function plotmeniscus(m,numfig,plotstyle);
            % function to plot the meniscus profile
            figure(numfig);
            plot(m.R,m.Z,plotstyle);
            hold on;
        end 
        
        
    end
end