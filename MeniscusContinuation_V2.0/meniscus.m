classdef meniscus
    %%
    %% Package to compute the equilibrium shapes and stability properties
    %% of hanging drops, attached bubbles, liquid bridges, and many other
    %% menisci problems controled by Young-Laplace equation
    %%
    %% Version 2.5 (D. Fabre, may 2015)
    %% Written in object-oriented matlab style
    %%
    
    % remarque : remplacer les (if (typestart == ''...)
    % par switch 
    
    properties
        %%%% physical properties
        nbdim = 3
        rhog = -1 % set to +1 for hanging drops and liquid bridges ; -1 for sessile drops
        gamma = 1 % warning : this should be introduced in the code in a clean way
        typestart = 'axisX' % possible types : 'axisX', 'pined' 
        typeend = 'pined'   % possible types : 'pined', 'angle'
        beta = 0            % contact angle, for angle-controled case
        Vref=10             % reference volume, for nondimensionalization and arclength;
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
        %%%% deformations in the iteration processes
        eta
        dP
        xi
        %%%% numerical parameters
        Vans
        Pans
        Vtarget
        discretization = 'FD'; %other choice 'FE'
        conv = 0;
        verbosity = 5;
        Weta = 1;
        Wxi = 1;
        typecont = 'P';
        % variable for continuation : 'P' (pressure), "V" (volume) or "S" (arclength)
        % Additional parameters for arclength continuation :
        dS;
        dVdS = 0;
        dPdS = 1;
        Nstep = 0;
        istab = 'no'; % if 'yes' we will compute and plot the stability properties
        ifig20 = 'no'; % if 'yes' plot additional results in figures 21,22,23,26
        idebug = 0;% set to larger values for 'debug' mode
        ixi = 0; % set to 1 to use grid reuniformisation using 'calcxi'
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% disp('diag centrale (A)');
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function m = meniscus(type,N,param)
            % Constructor for a guess Meniscus
            % type = flat => param is Radius
            % type = cylinder => param is [Radius,Length,Bond]
            % type = sphere => spherical cap, angle-controled, param is [R0,beta]
            
            
            switch(type)
            
                case('flat')
                m.N = N;
                R0 = param(1);
                t = 0:1/(N-1):1;
                m.R = param*t;
                m.Z = 0*t;
                m.V = 0;
                m.P = 0;
                m.typestart = 'axisX';
                m.typeend = 'pined';
                m.dVdS = 1; % to initiate arclength continuation in increasing-V
                m.dPdS = 0;
                
                case ('cylinder')
                R0 = param(1);
                L = param(2);
                m.rhog = param(3);
                m.N = N;
                t = 0:1/(N-1):1;
                m.R = R0*ones(size(t));
                m.Z = t*L;
                m.V = pi*R0^2*L;
                m.Vref = m.V;
                m.P = 1;
                m.typestart = 'pined';
                m.typeend = 'pined';
                
                case ('sphere')
                    % hanging portion of sphere
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
                m.rhog = -1; % -1 for bubble at ceiling ; +1 for hanging drop
                m. discretization = 'FE';
                m.idebug = 100;
                m.dPdS = -1; % to initiate arclength continuation in decreasing-P direction
             
                case ('Sessilesphere')
                    % hanging portion of sphere
                R0 = param(1);
                m.beta = param(2);
                m.N = N;
                t = 1:-1/(N-1):0;
                t = (m.beta)*t;
                m.R = R0*sin(t);
                m.Z = -R0*(cos(t(1))-cos(t));
                m.P = 2/R0;
                m.typeend = 'axisX';
                m.typestart = 'angle';%'pined';
                m.rhog = +1; % -1 for sessile drop ; +1 for hanging drop
                m. discretization = 'FE';
                m.beta = pi-m.beta;
                m.idebug = 100;
                m.dPdS = -1; % to initiate arclength continuation in decreasing-P direction
                case default
                disp(' Meniscus construction : invalid type !');
            end %switch
            
            m = calcgeom(m); 
            plotmeniscus(m,11,'k:');
            
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

            
            while ((RES > epsilon)&&(it<itmax)&&(RES<100))
                % Newton iteration 
                
                m = calcgeom(m);
                if(max(abs(m.ds))>1)
                    RES = 20000;
                end
                
                m = calcmat(m);
                
                if (m.typecont == 'P')
                    m = calcRHS(m);
                    
                    % debug sauvage
                    %etaattendu = 1e-3*linspace(-0.5968,-0.7010,50)
                    %RHSverif = m.AP(1:m.N,1:m.N)*etaattendu'
                    %m.RHS
                    %fin debug sauvage
                    
                    size(m.AP)
                    size(m.RHS)
                    
                    m.eta = (m.AP(1:m.N,1:m.N)\m.RHS(1:m.N)')';
                     
                    if(m.idebug>=10) 
                    disp(['debug : iteration ',num2str(it)]);
                    disp('RHS :');
                    m.RHS
                    disp('eta :');
                    m.eta
                    disp('3 diag centrales (A)');
                    [[diag(m.AP,-1);0],diag(m.AP),[0;diag(m.AP,1)]]
                     disp('valeurs propres(A)');
                    eig(m.AP)
                    plotmeniscus(m,11,'k:');
                    pause;
                    end
     
                    dP = 0;
                    RES = max(abs(m.RHS(1:m.N)));
                end
                
                if (m.typecont == 'V')
                    m = calcRHS(m);
                    etaXX = (m.AV(1:m.N+1,1:m.N+1)\m.RHS(1:m.N+1)')';
                    m.eta = etaXX(1:m.N);
                    dP = etaXX(m.N+1);
                    
                    
                   if(m.idebug>=10) 
                    disp(['debug : iteration ',num2str(it)]);
                    disp('RHS :');
                    m.RHS
                    disp('eta :');
                    m.eta
                    disp('dp : ');
                    dP
                    plotmeniscus(m,11,'k:');
                    pause;
                    end
                    
                    
                    RES = max(abs(m.RHS(1:m.N+1)));
                end
             
                
                m.R = m.R + m.eta.*sin(m.alpha);
                m.Z = m.Z - m.eta.*cos(m.alpha);
                m.P = m.P+dP;
                 if (m.typeend=='angle')
                    ZN = m.Z(N);
                    m.Z = m.Z-ZN;
                 elseif (m.typestart=='angle')
                    Z1 = m.Z(1);
                    m.Z = m.Z-Z1;
                else
                    ZN = 0;
                end
                % correction to ensure equal spacing of gridpoints 
              if(m.ixi==1)  
               
                m = calcxi(m);
                m.R = m.R  + m.Wxi*m.xi.*cos(m.alpha);
                m.Z = m.Z  + m.Wxi*m.xi.*sin(m.alpha);
              end
                
              if(m.idebug>=10)
                   m=calcgeom(m);
                    disp('alpha :');
                    m.alpha
                    disp('Ka :');
                    m.Ka
                    disp('Kb : ');
                    m.Kb  
              end
                
                it = it+1;
                             
            end %Newton iteration 
           
             if(RES>10||it==100)
                fprintf(' Divergence in Newton iteration ! \n')
                m = mans;
                m.conv = 0;
                fprintf(' Previous (P,V) : ( %f , %f ) \n\n',m.P,m.V);
                fprintf('number of steps : %i ',m.Nstep);
            else
                m.conv = 1;
                if (m.verbosity > 0)
                    fprintf(' New converged meniscus shape for (P,V) = ( %f , %f ); it = %i \n',m.P,m.V,it);
                end
                plotmeniscus(m,11,'r-');                
                
                if (strcmp(m.ifig20,'yes')==1)
                figure(21);hold on;plot(m.P,m.V,'r.');
                figure(26);hold on;plot(m.P,m.V/m.Vref,'r.');
                
                figure(22);hold on;plot(m.P,m.Z(m.N) - m.Z(1),'r.');
                figure(23);hold on;plot(m.P,m.R(m.N),'r.');
                end
                
                
             end
            m.Nstep = m.Nstep+1;
        end % function step
        
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
         function m = loop(m,type,dX,Nstep);
              % Loop over Newton steps, to draw P/V diagrams. 
              % Including eigenvalue computation
              % continuation type may be 'P', 'V', 'S'
            m = step(m,type,0);
            Vtab = [m.V];
            Ptab = [m.P];
            plotmeniscus(m,10,'k');
            itloop = 0;
            m.typecont=type;
            while (itloop<Nstep && m.conv==1)
                m = step(m,type,dX);    
                if(strcmp(m.istab,'yes')==1)
                    m = stab(m);
                end 
                itloop = itloop+1;
                if(m.conv==1)
                    Vtab = [Vtab m.V];
                    Ptab = [Ptab m.P];
                end
            end
            
            figure(20);
            hold on;
            plot(Ptab,Vtab,'k');
            figure(26);
            hold on;
            plot(Ptab,Vtab/m.Vref,'k');
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
           % m.ds(N) = m.ds(N-1); % not used but to avoid possible bugs
            
            % angle at center of sides (intermediate array)
            for i=1:N-1
                alphaS(i) =  atan2(m.Z(i+1)-m.Z(i),m.R(i+1)-m.R(i));
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
            
            
            if(m.typestart=='pined')
                m.alpha(1) = 2*m.alpha(2)-m.alpha(3);
                m.Ka(1) = 2*m.Ka(2)-m.Ka(3);
                m.Kb(1) = 2*m.Kb(2)-m.Kb(3);
            elseif(m.typestart=='angle')
                m.alpha(1) =2*m.alpha(2)-m.alpha(3);
                m.Ka(1) =2*m.Ka(2)-m.Ka(3);
                m.Kb(1) =2*m.Kb(2)-m.Kb(3);
            else % axis
                m.alpha(1) = 0;%0.;
                m.Ka(1) = 2*atan2(m.Z(2)-m.Z(1),m.R(2)-m.R(1) )/m.ds(1); 
                if (m.nbdim == 3)
                    m.Kb(1) = m.Ka(1);
                end
            end
            
            if (m.typeend=='pined')
                m.alpha(N) =2*m.alpha(N-1)-m.alpha(N-2);
                m.Ka(N) =2*m.Ka(N-1)-m.Ka(N-2);
                m.Kb(N) =2*m.Kb(N-1)-m.Kb(N-2);
            elseif (m.typeend=='axisX') 
                m.alpha(N) = pi;
                m.Ka(N) = 2*atan2(m.Z(N-1)-m.Z(N),m.R(N-1)-m.R(N) )/m.ds(N-1); 
                if (m.nbdim == 3)
                    m.Kb(N) = m.Ka(N);
                end
            else % angle
                %m.alpha(N) = m.beta; % ne marche pas en elements finis !
                % on essaie autre chose :S
                m.alpha(N) =2*m.alpha(N-1)-m.alpha(N-2);
                m.Ka(N) =2*m.Ka(N-1)-m.Ka(N-2);
                m.Kb(N) =2*m.Kb(N-1)-m.Kb(N-2);
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
            % Equation : -K1(eta)+N0z*eta*rhog = (K0 - (P+rhog*z))
            % Resolution as a matricial equation : A * eta = RHS
           
            A = 0;
            N = m.N;
            ds = m.ds;
            R = m.R;
            Z = m.Z;
            alpha = m.alpha;
            Ka = m.Ka;
            Kb = m.Kb;
            switch (m.discretization)
                
                case('FD') % finite difference
                    
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
                    % + N0z*eta*rhog%
                    A(i,i) = A(i,i) + m.rhog*cos(alpha(i));
                end;

                % case i=1
                if (strcmpi(m.typestart,'axisX')==1) %for bubble/drop : point 1 is on symmetry axis
                    if (m.nbdim==3)
                        A(1,1) = -4/(ds(1)^2);  %%% factor 4 because 1/R d eta / ds = d^2 eta / d s^2
                        A(1,2) = 4/(ds(1)^2);
                    else
                        A(1,1) = -2/(ds(1)^2);
                        A(1,2) = 2/(ds(1)^2);
                    end
                    A(1,1) = A(1,1)+(Ka(1)^2+Kb(1)^2)+m.rhog;
                elseif(strcmpi(m.typestart,'pined')==1) % for bridge : point 1 is pinned
                    A(1,1) = -1;
                    A(1,2) = 0;
                elseif(strcmpi(m.typestart,'angle')==1)
                    A(1,1) =  1/(2*ds(1));
                    A(1,2) = -4/(2*ds(1));
                    A(1,3)   =  3/(2*ds(1)); 
                 else 
                    error('Error : typestart not recognized');
                end

                if (strcmpi(m.typeend,'pined')==1)
                    % "pinned" case
                    A(N,N) = -1;
                    A(N,N-1) = 0;
                    A(N,N-2) = 0;
                elseif(strcmpi(m.typeend,'angle')==1)
                    A(N,N-2) =  1/(2*ds(N-1));
                    A(N,N-1) = -4/(2*ds(N-1));
                    A(N,N)   =  3/(2*ds(N-1));
                      % case i=1
                elseif (strcmpi(m.typeend,'axisX')==1) %for bubble/drop : point 1 is on symmetry axis
                    if (m.nbdim==3)
                        A(N,N) = -4/(ds(N)^2);  %%% factor 4 because 1/R d eta / ds = d^2 eta / d s^2
                        A(N,N-1) = 4/(ds(N)^2);
                    else
                        A(N,N) = -2/(ds(N)^2);
                        A(N,N-1) = 2/(ds(N)^2);
                    end
                    A(N,N) = A(N,N)+(Ka(N)^2+Kb(N)^2)+m.rhog;
                else 
                    error('Error : typeend not recognized');
                end
                
                AA1 = A(1:N,1:N); % for matrix 1, before adding shift for angle-controled case



                 if(strcmpi(m.typeend,'angle')==1)
                     % this term is because we displace the whole interface
                     % by eta n - eta(smax) cos(alphamax)
                    for i=1:N-1
                        A(i,N) = A(i,N) - m.rhog*cos(alpha(N));
                    end
                end

                m.AP = A;
                m.AV(1:N,1:N) = A;


                %%% END OF MATRIX AP ; NEXT IS FOR AV 

                % last column  for AV  / AS : effet of dP
                for i=2:N-1
                    m.AV(i,N+1) = +1;
                end;
                if(R(1)==0)
                    m.AV(1,N+1) = +1;
                end
                m.AV(N+1,N+1) = 0;

                % last line for AV : Volume increment
                for i=2:N-1
                    m.AV(N+1,i) =  -(ds(i-1)/2*(2*R(i)+R(i-1))/3 + ds(i)/2*(2*R(i)+R(i+1))/3)*2*pi;
                end;
                m.AV(N+1,1) = -ds(1)/2*(2*R(1)+R(2))/3*2*pi;
                m.AV(N+1,N) = -ds(N-1)/2*(2*R(N)+R(N-1))/3*2*pi;

                if(strcmpi(m.typeend,'angle')==1)
                    m.AV(N+1,N) = m.AV(N+1,N)-pi*R(N)^2*cos(m.alpha(N));
                end

                % Matrix A1 : for non-axisymetric instabilities

                for i=2:N
                    AA1(i,i) = AA1(i,i) - 1/R(i)^2;
                end

                AA1(1,1) = -1; AA1(1,2) = 0;  
                % point i=1 is either pinned or axis => cannot move
                % NB. point i=N (pined or angle) => same as in AP, nothing to do
                m.A1 = AA1;
             %fin Calcul de A par FD
            
            %Methode des elements finis
                case('FE')
                    
                dpdz=m.rhog;
               for i=2:N-1 % pour construire la ligne i etabar = hat(eta)(i) 
                    %  int ( - (d etabar / d s0 * d eta /d s0) *R ) ds0
                    A(i,i-1) = (R(i)+R(i-1))/2*1/ds(i-1);
                    A(i,i) = - R(i)*(1/ds(i)+1/ds(i-1));% 
                    % NB CA SERAIT PLUS LOGIQUE D'AVOIR :
                    % A(i,i) = - ( 1/ds(i)*(R(i)+R(i+1))/2 + 1/ds(i-1)*(R(i)+R(i-1))/2 )  
                    %if (i~=N-1)&&(m.typeend=='pined') 
                    A(i,i+1) =  (R(i)+R(i+1))/2*1/ds(i)  ; 
                    %end
                    %  int (Ka^2 + Kb^2- dpdz*cos(alpha)) etabar eta 
                    A(i,i-1) = A(i,i-1)  + (  (Ka(i)^2+Kb(i)^2-dpdz*cos(alpha(i)))*R(i) ...
                                            + (Ka(i-1)^2+Kb(i-1)^2-dpdz*cos(alpha(i-1)))*R(i-1) )/2 * (ds(i-1)/6);  
                    A(i,i) = A(i,i) + (Ka(i)^2+Kb(i)^2- dpdz*cos(alpha(i)))*R(i)*(ds(i-1)+ds(i))/3; 
                    %if (i~=N-1)&&(m.typeend=='pined') 
                    A(i,i+1) = A(i,i+1)  + (  (Ka(i)^2+Kb(i)^2-dpdz*cos(alpha(i)))*R(i) ...
                                            + (Ka(i+1)^2+Kb(i+1)^2-dpdz*cos(alpha(i+1)))*R(i+1) )/2 * (ds(i)/6); 
                    %end
                    %if(m.typeend=='pined')
                    %    A(N,N-1) = 0; % not even sure it is necessary
                    %end
               end
               
               
               % case i=1
                switch(m.typestart)
                
                    case('axisX') %for bubble/drop : point 1 is on symmetry axis
                  A(1,1) = -1/ds(1)*(R(2)/2);  %%% facteur 4 car 1/R d eta / ds = d^2 eta / d s^2 (mais en fait 2...) 
                  A(1,2) = A(2,1);
                  A(1,1) = A(1,1)+ (   (Ka(1)^2+Kb(1)^2-dpdz*cos(alpha(1)))*R(1) ...
                                     + (Ka(2)^2+Kb(2)^2-dpdz*cos(alpha(2)))*R(2) )/2 * (ds(1)/3); 
                    case('pined') % for bridge : point 1 is pinned
                  A(1,1) = 1;
                  A(1,2) = 0;
                
                    case('angle')
                        % a completer et tester
                    A(1,2) = (R(1)+R(2))/2*1/ds(1);
                    A(1,1) = - R(1)*(1/ds(1));% A(N,1) = 0;% a verifier 
                    A(1,2) = A(1,2)  + (  (Ka(1)^2+Kb(1)^2-dpdz*cos(alpha(1)))*R(1) ...
                                            + (Ka(2)^2+Kb(2)^2-dpdz*cos(alpha(2)))*R(2) )/2 * (ds(1)/6);  
                    A(1,1) = A(1,1) + (Ka(1)^2+Kb(1)^2- dpdz*cos(alpha(1)))*R(1)*(ds(1))/3;
                    A(1,3) = 0;
                    
                    A(N,1) = 0;
                    % terme supplementaire traduisant la translation
                    % verticale pour recoller à la surface
                    for i=2:N-1
                    A(i,1) = A(i,1) - dpdz*cos(alpha(1))*R(i)*(ds(i-1)+ds(i))/2;
                    end; 
                    A(1,1) = A(1,1) - dpdz*cos(alpha(1))*(R(1)+R(2))/2*(ds(1))/2;
                    A(N,1) = A(N,1) - dpdz*cos(alpha(1))*(R(N)+R(N-1))/2*(ds(N-1))/2;
                end
                

                
                
                switch(m.typeend)
                    
                    case('pined')
                    % "pinned" case
                    A(N,N) = -1;
                    A(N,N-1) = 0;
                    A(N,N-2) = 0;
                    
                    case('angle')
                    % cas a traiter différemment  
                    % partie integration par partie idem au cas i
                    
                    A(N,N-1) = (R(N)+R(N-1))/2*1/ds(N-1);
                    A(N,N) = - R(N)*(1/ds(N-1));% 
                    A(N,N-1) = A(N,N-1)  + (  (Ka(N)^2+Kb(N)^2-dpdz*cos(alpha(N)))*R(N) ...
                                            + (Ka(N-1)^2+Kb(N-1)^2-dpdz*cos(alpha(N-1)))*R(N-1) )/2 * (ds(N-1)/6);  
                    A(N,N) = A(N,N) + (Ka(N)^2+Kb(N)^2- dpdz*cos(alpha(N)))*R(N)*(ds(N-1))/3;
                    A(N,N-2) = 0;%pour etre sur
                    % terme supplementaire venant du terme de bord de l'IPP
                    % est a mettre dans le RHS
                    
                    % terme supplementaire traduisant la translation
                    % verticale pour recoller à la surface
                    for i=2:N-1
                    A(i,N) = A(i,N) - dpdz*cos(alpha(N))*R(i)*(ds(i-1)+ds(i))/2;
                    end; 
                    A(1,N) = A(1,N) - dpdz*cos(alpha(N))*(R(1)+R(2))/2*(ds(1))/2;
                    A(N,N) = A(N,N) - dpdz*cos(alpha(N))*(R(N)+R(N-1))/2*(ds(N-1))/2;
                 
                    case('axisX') %for bubble/drop : point 1 is on symmetry axis
                  A(N,N) = -1/ds(N-1)*(R(N-1)/2);  %%% facteur 4 car 1/R d eta / ds = d^2 eta / d s^2 (mais en fait 2...) 
                  A(N,N-1) = A(N-1,N);
                  A(N,N) = A(N,N)+ (   (Ka(N)^2+Kb(N)^2-dpdz*cos(alpha(N)))*R(N) ...
                                     + (Ka(N-1)^2+Kb(N-1)^2-dpdz*cos(alpha(N-1)))*R(N-1) )/2 * (ds(N-1)/3); 
                    
                end
                
                
                
                m.AP = A;
                m.AV(1:N,1:N) = A;
                
                % derniere ligne : 1/(2*pi) dV 
                for i=2:N-1
                    m.AV(N+1,i) =  R(i)*(ds(i)+ds(i-1))/2; % test sign
                    % ca semble + logique d'avoir
                    % (R(i)+R(i-1)/2)*ds(i-1)/2 + (R(i)+R(i+1))/2*ds(i)/2  
                    % a tester... 
                end;
                
                if (strcmpi(m.typestart,'axisX')==1)
                    m.AV(N+1,1) = R(2)/2*ds(1)/2; % test sign
                end
                
                if(strcmpi(m.typeend,'pined')==1)
                       m.AV(N+1,N) = R(N-1)/2*ds(N-1)/2; 
                elseif (strcmpi(m.typeend,'axisX')==1)
                    m.AV(N+1,N) = 0*R(N-1)/2*ds(N-1)/2; % test sign
                else
                    m.AV(N+1,N) = 0*R(N-1)/2*ds(N-1)/2;
                end
                m.AV(N+1,N+1) = 0;

                % last column : effet of dP
                for i=1:N-1
                    m.AV(i,N+1) = m.AV(N+1,i);
                end

            end%switch
        
        end % function calcmat
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function m = calcRHS(m);
            % computes the right-hand-side in P,V and S continuations
            dpdz=m.rhog;
            % Laplace/Young operator for central points
            switch(m.discretization)
                
                case('FD')
                
                for i=2:m.N-1
                    m.RHS(i) = m.gamma*(m.Ka(i)+m.Kb(i))-m.P+m.rhog*(m.Z(i)-m.Z(m.N));
                end
                if (strcmpi(m.typestart,'axisX')==1)
                    m.RHS(1) = m.gamma*(m.Ka(1)+m.Kb(1))-m.P+m.rhog*(m.Z(1)-m.Z(m.N));
                end
                if (strcmpi(m.typestart,'pined')==1)
                    m.RHS(1) = 0;
                end

                if (strcmpi(m.typeend,'pined')==1)
                    m.RHS(m.N) = 0;
                end
                if (strcmpi(m.typeend,'angle')==1)
                    m.RHS(m.N)   =  (m.alpha(m.N)-m.beta);
                end

                % For V continuation : volume constraint
                if (m.typecont=='V')
                    m.RHS(m.N+1) = -(m.Vtarget-m.V);
                end   
                
                
                
                case('FE') %FE
                    
                for i=2:m.N-1
                    m.RHS(i) = (m.gamma*(m.Ka(i)+m.Kb(i))-m.P+dpdz*m.Z(i))*m.R(i)*(m.ds(i-1)+m.ds(i))/2; 
                end
                if (strcmpi(m.typestart,'axisX')==1)
                    m.RHS(1) =0;
                end
                if (strcmpi(m.typestart,'pined')==1)
                    m.RHS(1) = 0;
                end
                if (strcmpi(m.typestart,'angle')==1)
                    % a verifier
                    m.RHS(1) =-(m.alpha(1)-m.beta)*(m.R(1));
                end
                
                
                
                if (strcmpi(m.typeend,'pined')==1)
                    m.RHS(m.N) = 0;
                end
                if (strcmpi(m.typeend,'angle')==1)
                    % a verifier
                    m.RHS(m.N) =-(m.alpha(m.N)-m.beta)*(m.R(m.N));
                end
                 if (strcmpi(m.typeend,'axisX')==1)
                    m.RHS(m.N) =0;
                end

                % For V continuation : volume constraint (warning in FE
                % we actually impose V/2pi) 
                if (m.typecont=='V')
                    m.RHS(m.N+1) = (m.Vtarget-m.V)/(2*pi);
                end
                
            end%switch  
            


        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function m = calcGuess(m);
            % Guess for S-iteration using NewtonEuler method
            
                RHS  = -m.AV(1:m.N,m.N+1);
                m.eta = (m.AP\RHS)';
                dP = 1;
                dV = -m.AV(m.N+1,1:m.N)*m.eta'/m.Vref;
                dS = sqrt(dP^2+dV^2);
                m.R = m.R + m.dS*m.dPdS*m.eta.*sin(m.alpha);
                m.Z = m.Z - m.dS*m.dPdS*m.eta.*cos(m.alpha);
                m.P = m.P+m.dS*m.dPdS;
                m = calcgeom(m);       
                
        end % function calcGuess
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function m = calcxi(m);
            m = calcgeom(m);
            N = m.N;
            curvelength = m.s0(N);
            for i=1:N
                m.xi(i)=(i-1)/(N-1)*(curvelength)-m.s0(i);
            end
        end % function calcxi
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                  
        function m = stab(m);
            % computes the minimun eigenvalues of AP,AV,A1 and checks if a
            % bifurcation is encountered in the loop process
            
            lambdaminVans = m.lambdaminV;
            lambdaminPans = m.lambdaminP;
            lambdamin1ans = m.lambdamin1;
            
            EE = eig(m.AP);
            EE = EE(imag(EE)==0);
            EV = EE(abs(EE)==min(abs(EE)));
            m.lambdaminP = EV(1);
            
            EE = eig(m.AV);
            EE = EE(imag(EE)==0);
            EV = EE(abs(EE)==min(abs(EE)));
            m.lambdaminV = EV(1); 
            
            EE = eig(m.A1);
            EE = EE(imag(EE)==0);
            EV = EE(abs(EE)==min(abs(EE)));
            m.lambdamin1 = EV(1);
            
            if(m.verbosity>1)
                figure(30);
                hold on;
                plot(m.P,real(m.lambdaminP),'r+',m.P,real(m.lambdaminV),'b+',m.P,real(m.lambdamin1),'g+');
                plot(m.P,imag(m.lambdaminP),'rx',m.P,imag(m.lambdaminV),'bx',m.P,imag(m.lambdamin1),'gx');
            end
            
            % detection bifurcation P
            if (real(m.lambdaminP)*real(lambdaminPans)<0&&abs(m.lambdaminP-lambdaminPans)<1&&abs(imag(m.lambdaminP))<0.5) 
                VbV = (m.V*lambdaminPans - m.Vans*m.lambdaminP)/(lambdaminPans-m.lambdaminP);
                PbV = (m.P*lambdaminPans - m.Pans*m.lambdaminP)/(lambdaminPans-m.lambdaminP);
                figure(20);
                hold on;
                plot(PbV,VbV,'xr','MarkerSize',10);
                figure(25);
                hold on;
                plot(PbV,VbV/m.Vref,'xr','MarkerSize',10);
                if (imag(VbV)==0)
                    mbif = step(m,'V',VbV-m.V);
                    plotmeniscus(mbif,10,'r');
                end
                fprintf('***\n*** P-bifurcation detected for (P,V) = ( %f , %f ) \n***\n', PbV,VbV);
            end
            
            % detection bifurcation V
             if (real(m.lambdaminV)*real(lambdaminVans)<0&&abs(m.lambdaminV-lambdaminVans)<1&&abs(imag(m.lambdaminV))<0.5) 
                VbV = (m.V*lambdaminVans - m.Vans*m.lambdaminV)/(lambdaminVans-m.lambdaminV);
                PbV = (m.P*lambdaminVans - m.Pans*m.lambdaminV)/(lambdaminVans-m.lambdaminV);
                figure(20);
                hold on;
                plot(real(PbV),real(VbV),'+b','MarkerSize',10);
                figure(25);
                hold on;
                plot(real(PbV),real(VbV/m.Vref),'+b','MarkerSize',10);
                if (imag(VbV)==0)
                    mbif = step(m,'P',real(PbV)-m.P);
                    plotmeniscus(mbif,10,'b');
                 end
                 fprintf('***\n*** V-bifurcation detected for (P,V) = ( %f , %f ) \n***\n', PbV,VbV);
             end
            
            % detection bifurcation m=1
            if (m.lambdamin1*lambdamin1ans<0&&abs(m.lambdamin1-lambdamin1ans)<1&&imag(m.lambdamin1)==0) 
                VbV = (m.V*lambdamin1ans - m.Vans*m.lambdamin1)/(lambdamin1ans-m.lambdamin1);
                PbV = (m.P*lambdamin1ans - m.Pans*m.lambdamin1)/(lambdamin1ans-m.lambdamin1);
                figure(20);
                hold on;
                plot(PbV,VbV,'*g','MarkerSize',10);
                figure(25);
                hold on;
                plot(PbV,VbV/m.Vref,'*g','MarkerSize',10);
                if (imag(VbV)==0)
                    mbif = step(m,'V',VbV-m.V);
                    plotmeniscus(mbif,10,'g');
                end
                fprintf('***\n*** NonAxi-bifurcation detected for (P,V) = ( %f , %f ) \n***\n', PbV,VbV);
            end
            
          
            
        end % function stab
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        function plotmeniscus(m,numfig,plotstyle);
            % function to plot the meniscus profile
            figure(numfig);
            plot(m.R,m.Z,plotstyle);
            hold on;
            grid on
        end %function plotmeniscus
        
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function m = flip(m);
            % to flip a meniscus profile (useful for bridges)
            mans = m;
            L = -m.Z(m.N);
            for i=1:m.N
                m.R(i) = mans.R(m.N+1-i);
                m.Z(i) = -L-mans.Z(m.N+1-i);
            end
            m = step(m,'P',0);
        end % function flip
        
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function resetfigs(m);
            % to reset figures
 
            close all;
            figure(10);
            title('Menisci shapes');
            xlabel('R');
            ylabel('Z');
            axis equal;
            
            hold on;
            
            figure(11);
            title('Menisci shapes (all computed)');
            xlabel('R');
            ylabel('Z');
            hold on;
            
            figure(20);
            title('P/V curves');
            xlabel('P');
            ylabel('V');
            hold on;
            
            if (strcmp(m.ifig20,'yes')==1)
            
            figure(21);
            title('P/V results (dots)');
            xlabel('P');
            ylabel('V');
            hold on;
            
            figure(22);
            title('P/H results');
            xlabel('P');
            ylabel('H');
            hold on;
            
            figure(23);
            title('P/R0 results');
            xlabel('P');
            ylabel('R0');
            hold on;
            
            
            figure(26);
            title('P/V* results');
            xlabel('P');
            ylabel('V*');
            hold on;
            
            end
            
            if(strcmp(m.istab,'yes')==1)
            figure(30);
            title('Eigenvalues of matrices AP, AV, A1');
            ylabel('lambdaP, lambdaV, lambda1');
            xlabel('P');
            hold on;
            end
            
        end % function resetfigs
           
    end % methods

end % class meniscus 