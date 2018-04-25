classdef meniscus
    %%
    %% Package to compute the equilibrium shapes and stability properties
    %% of hanging drops, attached bubbles, liquid bridges, and many other
    %% menisci problems controled by Young-Laplace equation
    %%
    %% Version 2.5 (D. Fabre, may 2015)
    %% Written in object-oriented matlab style
    %%
    
   
    
    properties
        %%%% physical properties
        nbdim = 3
        rhog = 1  % set to 1 to work in nondimensional units
        gamma = 1 % idem
        typestart = 'axisX' % possible types : 'axisX', 'pined' 
        typeend = 'pined'   % possible types : 'pined', 'angle'
        beta = 0            % contact angle, OLD VERSION 
        alphastart          % if typestart = 'angle', this will be the imposed value of the angle; in other cases this parameter is not used.
        alphaend            % same thing if typeend = 'angle'
        Vref=3             % reference volume, for nondimensionalization and arclength;
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
        lambdaminV = -1;
        lambdaminP = -1;
        lambdamin1 = -1;
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
        dVdS = 1;
        dPdS = 0;
        Nstep = 0;
        istab = 'no'; % if 'yes' we will compute and plot the stability properties
        whichfigures = [10 20] % selection of figures to plots. (10 = shapes; 20 = PV diagram; many more available...)
        idebug = 0;% set to larger values for 'debug' mode
        ixi = 1; % set to 1 to use grid reuniformisation using 'calcxi'
        isinloop = 0; % this is an internal parameter to disable plots when doing loops
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
               
                case('flat') % use this to initiate a sessile drop
                m.N = N;
                R0 = param(1);
                t = 1:-1/(N-1):0;
                m.R = param*t;
                m.Z = 0*t;
                m.V = 0;
                m.P = 0;
                m.typestart = 'pined';
                m.typeend = 'axisX';
                m.rhog = +1;
                
                case('flatinv') % use this to initiate a hanging drop
                m.N = N;
                R0 = param(1);
                t = 0:1/(N-1):1;
                m.R = param*t;
                m.Z = 0*t;
                m.V = 0;
                m.P = 0;
                m.typestart = 'axisX';
                m.typeend = 'pined';
                m.rhog = +1;
                m.verbosity = 10;
                
                case ('cylinder') % use this to initiate a liquid bridge
                R0 = param(1);
                L = param(2);
                m.rhog = param(3);
                m.N = N;
                t = 0:1/(N-1):1;
                m.R = R0*ones(size(t));
                m.Z = t*L;
                m.V = pi*R0^2*L;
                m.Vref = m.V;
                m.P = m.gamma/R0;
                m.typestart = 'pined';
                m.typeend = 'pined';
                
                case ('cap') % spherical cap posed onto a plate
                m.N = N;
                
                if(param(2)<=pi)
                    error('WARNING : value of angle seems very small... is that the right value in DEGREES ????');
                end
                thetas = param(2)*pi/180; % static contact angle (converted to radians)
                R0 = param(1)/sin(thetas); % this is the radius AT THE CONTACT LINE 
                m.typeend = 'axisX';
                m.typestart = 'angle';%'pined';
                m.alphastart = pi-thetas;
                t = 1:-1/(N-1):0;
                t = thetas*t;
                m.R = R0*sin(t);
                m.Z = -R0*(cos(thetas)-cos(t));                
                m.P = 2*m.gamma/R0; % initial pressure if gravity is negliged
                m.rhog = +1; % 
                m. discretization = 'FE';
                m.beta = pi-thetas;
                m.dPdS = 1; % to initiate arclength continuation in increasing-P direction

                case ('capinv') % spherical cap hanging from a plate
                m.N = N;
                 if(param(2)<=pi)
                    error('WARNING : value of angle seems very small... is that the right value in DEGREES ????');
                end
                thetas = param(2)*pi/180; % static contact angle (converted to radians)
                 R0 = param(1)/sin(thetas); % this is the radius AT THE CONTACT LINE
                m.typestart = 'axisX';
                m.typeend = 'angle';
                m.alphaend = thetas;
                t = 0:1/(N-1):1;
                t = thetas*t;
                m.R = R0*sin(t);
                m.Z = R0*(cos(thetas)-cos(t));                
                m.P = 2*m.gamma/R0; % initial pressure if gravity is negliged
                m.rhog = +1; % 
                m. discretization = 'FE';
                m.beta = thetas;
                m.dPdS = -1; % to initiate arclength continuation in decreasing-P direction
    
                
                case ('sphere') % previous method for constructing a hanging portion of sphere
                R0 = param(1);
                m.beta = param(2);
                R0 = R0/sin(m.beta); % R0 is the initial radius 
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
             
                case ('Sessilesphere') % previous method for constructing a sessile portion of sphere
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
                error(' Meniscus construction : invalid type !');
            end %switch
            
            m.resetfigs; % to set the axes and legends for the plots
            m = calcgeom(m); 
            plotmeniscus(m,10,'k:');
            
        end % function meniscus
        
       
        
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        function m = step(mans,type,dX)
            % Newton step from a previous meniscus.
            % continuation type may be 'P', 'V', 'dP', 'dV', 'dS'
            % ('theta' and 'dtheta' should be implemented someday)
            
            m = mans; 
            m.Vans = m.V;
            m.Pans = m.P;
            
            if(strcmp(m.discretization,'FE')==1&&strcmp(m.istab,'yes')==1)
                error('ERROR in meniscus/step : stability not yet implemented under FE discretization (use FD)');
            end
            if(strcmp(m.discretization,'FE')==1&&strcmp(type,'dS')==1)
                error('ERROR in meniscus/step : arclength continuation not yet implemented under FE discretization (use FD)');
            end
            
           
            m = calcgeom(m);
            m.conv = 0; 
            N = m.N;
            
            RES = 1;
            epsilon = 1e-6;
            it = 0;
            itmax = 100;
  
            switch(type)    
            case('dP')
                m.P = m.P+dX;
                m.typecont = 'P';   
            case('dV')
                m.Vtarget = m.V+dX;
                m.typecont = 'V';
            case('P')
                m.P = dX;
                 m.typecont = 'P';
            case('V')
                m.Vtarget = dX;
                m.typecont = 'V';
            case('dS')
                m.typecont = 'S'; 
                m = calcmat(m);
                m.dS = dX;
                m = calcGuess(m);
                  if(m.idebug>=10) 
                    disp(['debug : guess pour S ',num2str(it)]);
                    m.R 
                    m.Z
                    plotmeniscus(m,11,'b:');
                    pause;
                  end
%            case default 
%                    error(' ERROR in meniscus/step : type of loop not recognized (use P,dP,V,dV or S)');
            end %switch 

            
            while ((RES > epsilon)&&(it<itmax)&&(RES<100))
                % Newton iteration 
                
                m = calcgeom(m);
                if(max(abs(m.ds))>1)
                    RES = 20000;
                end
                
                m = calcmat(m);
                
                switch(m.typecont)
                    case('P')
                    m = calcRHS(m);
                    m.eta = (m.AP(1:m.N,1:m.N)\m.RHS(1:m.N)')';                  
                    dP = 0;
                    RES = max(abs(m.RHS(1:m.N)))/max(m.gamma*max(m.R)+abs(m.rhog)*max(abs(m.Z))); 
                case('V')
                    m = calcRHS(m);
                    etaXX = (m.AV(1:m.N+1,1:m.N+1)\m.RHS(1:m.N+1)')';
                    m.eta = etaXX(1:m.N);
                    dP = etaXX(m.N+1);
                    RES = max(abs(m.RHS(1:m.N+1)))/max(m.gamma*max(m.R)+abs(m.rhog)*max(abs(m.Z)));  
                case('S')
                    m = calcRHS(m);
                    etaXX = (m.AS(1:m.N+1,1:m.N+1)\m.RHS(1:m.N+1)')';
                    m.eta = etaXX(1:m.N);
                    dP = etaXX(m.N+1);
                    RES = max(abs(m.RHS(1:m.N+1)))/max(m.gamma*max(m.R)+abs(m.rhog)*max(abs(m.Z)));
                end%switch

                if(m.idebug>=10) 
                    disp(['debug : iteration ',num2str(it)]);
                    disp('RHS :');
                    m.RHS
                    %disp('3 diag centrales (A)');
                    %[[diag(m.AP,-1);0],diag(m.AP),[0;diag(m.AP,1)]]
                    m.AP
                    disp('valeurs propres(A)');
                    eig(m.AP)
                    disp('eta :');
                    m.eta
                    plotmeniscus(m,11,'k:');
                    pause;
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
                
            
                
                it = it+1;
                             
            end %Newton iteration 
           
             if(RES>10||it==100)
                fprintf(' Divergence in Newton iteration ! \n')
                m = mans;
                m.conv = 0;
                fprintf(' Previous (P,V) : ( %f , %f ) \n\n',m.P,m.V);
                fprintf('number of steps : %i \n \n',m.Nstep);
             else
                    
                m.conv = it;
               
                if(strcmp(m.typestart,'pined')==1)
                    m.beta = pi-m.alpha(1);
                end
                if(strcmp(m.typeend,'pined')==1)
                    m.beta = m.alpha(end);
                end
                 if (m.verbosity > 0)
                    fprintf(' New converged meniscus shape for (P,V,theta) = ( %g , %g , %g); it = %i \n',m.P,m.V,m.beta*180/pi,it);
                end
                
                
                %%% plots in figures 10 and 20 if not in loop
                 if(m.isinloop==0) 
                    plotmeniscus(m,10,'k-'); % plot in figure 10 only when step is called directly, not inside loop
                     if(ismember(20,m.whichfigures)||m.idebug>49)   
                     figure(20);hold on;plot(m.P,m.V,'k.');
                     end
                 end
                  
                plotmeniscus(m,11,'r-');                
                
               %%% MISCLEANOUS PLOTS ON DEMAND
                if(ismember(21,m.whichfigures)||m.idebug>49)   
                figure(21);hold on;plot(m.P,m.V,'r.');
                end
                if(ismember(26,m.whichfigures)||m.idebug>49)
                figure(26);hold on;plot(m.P,m.V/m.Vref,'r.');
                end
                if(ismember(22,m.whichfigures)||m.idebug>49)
                figure(22);hold on;plot(m.P,m.Z(m.N) - m.Z(1),'r.');
                end
                if(ismember(23,m.whichfigures)||m.idebug>49)
                figure(23);hold on;plot(m.P,m.R(m.N),'r.');
                end
                
                if(ismember(28,m.whichfigures)||m.idebug>49)
                figure(28);hold on;plot(max(m.R(1),m.R(m.N)),abs(m.Z(m.N)-m.Z(1)),'r.');
                end
                
             end
            m.Nstep = m.Nstep+1;
        end % function step
        
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
         function m = loop(m,type,dX,Nstep);
              % Loop over Newton steps, to draw P/V diagrams. 
              % Including eigenvalue computation
              % continuation type may be 'P', 'V', 'S'
             if (m.verbosity > 0)
                    disp([' Starting loop of type ',type,' with parameters ',num2str(dX),' ; ',num2str(Nstep)]) ;
             end
            m.isinloop=1; % to disable plots of all menisci shapes in figure 10    
            m = step(m,type,0);
            Vtab = [m.V];
            Ptab = [m.P];
            thetatab = [m.beta];
            plotmeniscus(m,10,'k');
            itloop = 0;
            
  %          switch(type)
  %              case('P')
  %          m.typecont='dP';
  %              case('V')
  %          m.typecont='dV';
  %              case('S')
  %          m.typecont='dS';
%                case(default) 
%                    error(' ERROR in meniscus/loop : type of loop not recognized (use P,V or S)');
  %          end
  m.typecont = type;          
            
            while (itloop<Nstep && m.conv>0)
                m = step(m,type,dX);    
                if(strcmp(m.istab,'yes')==1)
                    m = stab(m);
                end 
                itloop = itloop+1;
                if(m.conv>1)
                    Vtab = [Vtab m.V];
                    Ptab = [Ptab m.P];
                    thetatab = [thetatab m.beta];
                end
            end
            
            if (m.verbosity > 0 && m.conv>1)
                    disp([' Loop correctly completed']) ;
            end 
            
            if(ismember(20,m.whichfigures)||m.idebug>49)
            figure(20);
            hold on;
            plot(Ptab,Vtab,'k');
            end
            if(ismember(26,m.whichfigures)||m.idebug>49)
            figure(26);
            hold on;
            plot(Ptab,Vtab/m.Vref,'k');
            end
            
            if(ismember(121,m.whichfigures)||m.idebug>49)
            figure(121);
            hold on;
            plot(Vtab,thetatab,'k');
            end
            
            m.conv=1;
            plotmeniscus(m,10,'k');
            m.isinloop=0; % to re-able plots of all menisci shapes in figure 10   
            
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
            alphaS = unwrap(alphaS);
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
            
            
            switch(m.typestart)
                case('pined')
%                m.alpha(1) = 2*m.alpha(2)-m.alpha(3);
%                m.Ka(1) = 2*m.Ka(2)-m.Ka(3);
%                m.Kb(1) = 2*m.Kb(2)-m.Kb(3);
                alphaA = alphaS(1);alphaB = alphaS(2);alphaC = alphaS(3);
                sA = m.ds(1)/2; sB = m.ds(1)+m.ds(2)/2; sC = m.ds(1)+m.ds(2)+m.ds(3)/2;
                m.alpha(1) = (alphaA*sB*sB*sC-alphaA*sB*sC*sC-alphaB*sA*sA*sC+alphaB*sA*sC*sC+...
                             alphaC*sA*sA*sB-alphaC*sA*sB*sB)/(sA-sB)/(sA*sB-sA*sC-sB*sC+sC*sC);
                m.Ka(1) = -(alphaA*sB*sB-sC*sC*alphaA-sA*sA*alphaB+sC*sC*alphaB+sA*sA*alphaC-...
                            alphaC*sB*sB)/(sA-sB)/(sA*sB-sA*sC-sB*sC+sC*sC);
                m.Kb(1) = sin(m.alpha(1))/m.R(1);        
                case ('angle')
               % m.alpha(1) =2*m.alpha(2)-m.alpha(3);
               % m.Ka(1) =2*m.Ka(2)-m.Ka(3);
               % m.Kb(1) =2*m.Kb(2)-m.Kb(3);
                alphaA = alphaS(1);alphaB = alphaS(2);alphaC = alphaS(3);
                sA = m.ds(1)/2; sB = m.ds(1)+m.ds(2)/2; sC = m.ds(1)+m.ds(2)+m.ds(3)/2;
                m.alpha(1) = (alphaA*sB*sB*sC-alphaA*sB*sC*sC-alphaB*sA*sA*sC+alphaB*sA*sC*sC+...
                             alphaC*sA*sA*sB-alphaC*sA*sB*sB)/(sA-sB)/(sA*sB-sA*sC-sB*sC+sC*sC);
                m.Ka(1) = -(alphaA*sB*sB-sC*sC*alphaA-sA*sA*alphaB+sC*sC*alphaB+sA*sA*alphaC-...
                            alphaC*sB*sB)/(sA-sB)/(sA*sB-sA*sC-sB*sC+sC*sC);
                m.Kb(1) = sin(m.alpha(1))/m.R(1); 
                case('axisX')
                m.alpha(1) = 0;
                %m.Ka(1) = 2*atan2(m.Z(2)-m.Z(1),m.R(2)-m.R(1) )/m.ds(1);
                m.Ka(1) = 2*alphaS(1)/(m.ds(1));
                if (m.nbdim == 3)
                    m.Kb(1) = m.Ka(1);
                end
            case default
                    error('Error : typestart not recognized');
            end    
            
            switch(m.typeend)
                case('pined')
                m.alpha(N) =2*m.alpha(N-1)-m.alpha(N-2);
                m.Ka(N) =2*m.Ka(N-1)-m.Ka(N-2);
                m.Kb(N) =2*m.Kb(N-1)-m.Kb(N-2);
                case('axisX') 
                m.alpha(N) = pi; % actually \pm pi, the right sign is selected
                                 % by the unwrap command done below. 
                %m.Ka(N) = -2*atan2(m.Z(N-1)-m.Z(N),m.R(N-1)-m.R(N) )/m.ds(N-1);
                m.Ka(N) = 2*sin(alphaS(N-1))/m.ds(N-1); 
                % this should actually be (alphaS(N)-alphaS(N-1))/ds(N-1)
                % but alphaS(N) does not exist, and should be symmetric of
                % alphaS(N-1) \pm pi, and the sign may vary here !
                if (m.nbdim == 3)
                    m.Kb(N) = m.Ka(N);
                end
                case('angle')
                %m.alpha(N) = m.beta; % ne marche pas en elements finis !
                % on essaie autre chose :S
                m.alpha(N) =2*m.alpha(N-1)-m.alpha(N-2);
                m.Ka(N) =2*m.Ka(N-1)-m.Ka(N-2);
                m.Kb(N) =2*m.Kb(N-1)-m.Kb(N-2);
                case default
                    error('Error : typeend not recognized');
            end
            
            m.alpha = unwrap(m.alpha);
            
            m.V = 0;
            for i=2:N
                m.V = m.V+ ( (m.R(i)+m.R(i-1))^2/4 + (m.R(i)-m.R(i-1))^2/12 )*(m.Z(i)-m.Z(i-1))*pi;
            end;
            
              if(m.idebug>=10)
                   disp('GEOMETRIE : ');
                    disp('alpha :');
                    m.alpha
                    disp('Ka :');
                    m.Ka
                    disp('Kb : ');
                    m.Kb  
                    m.idebug
                    pause;
              end
            
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
                    A(i,i-1) = m.gamma*(2*ds(i)/(ds(i)*ds(i-1)*(ds(i)+ds(i-1))));
                    A(i,i) = m.gamma*(-2/(ds(i)*ds(i-1)));
                    A(i,i+1) = m.gamma*(2*ds(i-1)/(m.ds(i)*ds(i-1)*(ds(i)+ds(i-1))));
                    % T0r/r* d eta / d s0
                    if (m.nbdim==3)
                        A(i,i-1) = A(i,i-1) + m.gamma*(-ds(i)^2/((ds(i)+ds(i-1))*ds(i)*ds(i-1))*cos(alpha(i))/R(i));
                        A(i,i) = A(i,i) +  m.gamma*((ds(i)^2-ds(i-1)^2)/((ds(i)+ds(i-1))*ds(i)*ds(i-1))*cos(alpha(i))/R(i));
                        A(i,i+1) = A(i,i+1) +m.gamma*(ds(i-1)^2/((ds(i)+ds(i-1))*ds(i)*ds(i-1))*cos(alpha(i))/R(i));
                    end
                    %  (Ka^2 + Kb^2) eta
                    A(i,i) = A(i,i) + m.gamma*(Ka(i)^2+Kb(i)^2);
                    % + N0z*eta*rhog%
                    A(i,i) = A(i,i) + m.rhog*cos(alpha(i));
                end;

        % first point
                switch(m.typestart)
                    case('axisX') 
                    if (m.nbdim==3)
                        A(1,1) = m.gamma*(-4/(ds(1)^2));  %%% factor 4 because 1/R d eta / ds = d^2 eta / d s^2
                        A(1,2) = m.gamma*(4/(ds(1)^2));
                    else
                        A(1,1) = m.gamma*(-2/(ds(1)^2));
                        A(1,2) = m.gamma*(2/(ds(1)^2));
                    end
                    A(1,1) = A(1,1)+m.gamma*(Ka(1)^2+Kb(1)^2)+m.rhog*cos(alpha(1));
                    case('pined') 
                    A(1,1) = -1;
                    A(1,2) = 0;
                    case('angle')
                    A(1,1) =  1/(2*ds(1));
                    A(1,2) = -4/(2*ds(1));
                    A(1,3)   =  3/(2*ds(1)); 
                end
                
                % last point
                switch(m.typeend)
                    case('pined')
                    % "pinned" case
                    A(N,N) = -1;
                    A(N,N-1) = 0;
                    A(N,N-2) = 0;
                    case('angle')
                    A(N,N-2) =  1/(2*ds(N-1));
                    A(N,N-1) = -4/(2*ds(N-1));
                    A(N,N)   =  3/(2*ds(N-1));
                    case('axisX') %for bubble/drop : point 1 is on symmetry axis
                    if (m.nbdim==3)
                        A(N,N) = m.gamma*(-4/(ds(N-1)^2));  %%% factor 4 because 1/R d eta / ds = d^2 eta / d s^2
                        A(N,N-1) = m.gamma*(4/(ds(N-1)^2));
                    else
                        A(N,N) =m.gamma*( -2/(ds(N-1)^2));
                        A(N,N-1) = m.gamma*(2/(ds(N-1)^2));
                    end
                    A(N,m.N) = A(N,N)+m.gamma*(Ka(N)^2+Kb(N)^2)+m.rhog*cos(alpha(N));
                end
                
                AA1 = A(1:N,1:N); % for matrix 1, before adding shift for angle-controled case

                 if(strcmpi(m.typeend,'angle')==1)
                     % this term is because we displace the whole interface
                     % by eta n - eta(smax) cos(alphamax)
                    for i=1:N
                        A(i,N) = A(i,N) - m.rhog*cos(alpha(N));
                    end
                end

                m.AP = A;
                m.AV(1:N,1:N) = A;
                m.AS(1:N,1:N) = A;

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

                % Matrix AS : increment in arclength
            
                    m.AS(1:N,N+1) = m.AV(1:N,N+1);
                    m.AS(N+1,1:N) = -m.dVdS*m.AV(N+1,1:N)/m.Vref; 
                    m.AS(N+1,N+1) = m.dPdS;
                
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
                    %  int ( - gamma*(d etabar / d s0 * d eta /d s0) *R ) ds0
                    A(i,i) = - m.gamma*( 1/ds(i)*(R(i)+R(i+1))/2 + 1/ds(i-1)*(R(i)+R(i-1))/2 ); 
                    A(i,i-1) = m.gamma*(R(i)+R(i-1))/2*1/ds(i-1);
                    A(i,i+1) =  m.gamma*((R(i)+R(i+1))/2*1/ds(i))  ; 
             
                    %  int (Ka^2 + Kb^2- dpdz*cos(alpha)) etabar eta 
                   
                  %  A(i,i-1) = A(i,i-1)  + (  (m.gamma*(Ka(i)^2+Kb(i)^2)-dpdz*cos(alpha(i)))*R(i) ...
                  %                          + (m.gamma*(Ka(i-1)^2+Kb(i-1)^2)-dpdz*cos(alpha(i-1)))*R(i-1) )/2 * (ds(i-1)/6);   
                  %  A(i,i) = A(i,i) +  (  (m.gamma*(Ka(i)^2+Kb(i)^2)-dpdz*cos(alpha(i)))*R(i) ...
                  %                          + (m.gamma*(Ka(i-1)^2+Kb(i-1)^2)-dpdz*cos(alpha(i-1)))*R(i-1) )/2 * (ds(i-1)/3)... 
                  %  + (  (m.gamma*(Ka(i)^2+Kb(i)^2)-dpdz*cos(alpha(i)))*R(i) ...
                  %                          + (m.gamma*(Ka(i+1)^2+Kb(i+1)^2)-dpdz*cos(alpha(i+1)))*R(i+1) )/2 * (ds(i)/3); 
                  %  A(i,i+1) = A(i,i+1)  + (  (m.gamma*(Ka(i)^2+Kb(i)^2)-dpdz*cos(alpha(i)))*R(i) ...
                  %                         + (m.gamma*(Ka(i+1)^2+Kb(i+1)^2)-dpdz*cos(alpha(i+1)))*R(i+1) )/2 * (ds(i)/6);
                    KKi = m.gamma*(Ka(i)^2+Kb(i)^2)-dpdz*cos(alpha(i));  
                    dpKKi =  m.gamma*(Ka(i+1)^2+Kb(i+1)^2)-dpdz*cos(alpha(i+1))-KKi;
                    dmKKi =  m.gamma*(Ka(i-1)^2+Kb(i-1)^2)-dpdz*cos(alpha(i-1))-KKi; 
                    ri = R(i) ; dpri = R(i+1)-ri ; dmri = R(i-1)-ri ;
                    A(i,i) = A(i,i) + (20*KKi*ri+5*KKi*dpri+5*dpKKi*ri+2*dpKKi*dpri)*ds(i)/60 ...
                                    + (20*KKi*ri+5*KKi*dmri+5*dmKKi*ri+2*dmKKi*dmri)*ds(i-1)/60;
                    A(i,i-1) = A(i,i-1) + (10*KKi*ri+5*KKi*dmri+5*dmKKi*ri+3*dmKKi*dmri)*ds(i-1)/60;
                    A(i,i+1) = A(i,i+1) + (10*KKi*ri+5*KKi*dpri+5*dpKKi*ri+3*dpKKi*dpri)*ds(i)/60;
               end
               
               %A(1,1) = m.gamma*(-1/ds(1)*(R(1)+R(2))/2);
               KKi = m.gamma*(Ka(1)^2+Kb(1)^2)-dpdz*cos(alpha(1));  
               dpKKi =  m.gamma*(Ka(2)^2+Kb(2)^2)-dpdz*cos(alpha(2))-KKi;
               ri = R(1) ; dpri = R(2)-ri; 
               A(1,1) = A(1,1) + (20*KKi*ri+5*KKi*dpri+5*dpKKi*ri+2*dpKKi*dpri)*ds(1)/60;
               A(1,2) = A(2,1);
               
                switch(m.typestart)
                
                    case('axisX') %for hanging drop : point 1 is on symmetry axis 
                        % (nothing to add)
                 % A(1,1) = m.gamma*(-1/ds(1)*(R(1)+R(2))/2);  
                 % A(1,2) = A(2,1);
 %                 A(1,1) = A(1,1)+ (   (m.gamma*(Ka(1)^2+Kb(1)^2)-dpdz*cos(alpha(1)))*R(1) ...
 %                                    + (m.gamma*(Ka(2)^2+Kb(2)^2)-dpdz*cos(alpha(2)))*R(2) )/2 * (ds(1)/3); 
                 % A(1,1) = A(1,1)+ (   (m.gamma*(Ka(1)^2+Kb(1)^2)-dpdz*cos(alpha(1))) ...
                  %                   + (m.gamma*(Ka(2)^2+Kb(2)^2)-dpdz*cos(alpha(2))))*(R(1)+R(2))/2 * (ds(1)/3); 
                    case('pined') % for bridge : point 1 is pinned
                  A(1,1) = -1;
                  A(1,2) = 0;
                
                    case('angle') % for sessile drop with imposed angle
                   % A(1,2) = (R(1)+R(2))/2*1/ds(1);
                   % A(1,1) = - R(1)*(1/ds(1));% A(N,1) = 0;% a verifier 
                   % A(1,2) = A(1,2)  + (  (m.gamma*(Ka(1)^2+Kb(1)^2)-dpdz*cos(alpha(1)))*R(1) ...
                   %                         + (m.gamma*(Ka(2)^2+Kb(2)^2)-dpdz*cos(alpha(2)))*R(2) )/2 * (ds(1)/6);  
                   % A(1,1) = A(1,1) + (m.gamma*(Ka(1)^2+Kb(1)^2)- dpdz*cos(alpha(1)))*R(1)*(ds(1))/3;
                  
                    % terme supplementaire traduisant la translation
                    % verticale pour recoller ?? la surface
                    for i=2:N-1
                    A(i,1) = A(i,1) - dpdz*cos(alpha(1))*( (R(i)/2+(R(i+1)-R(i))/6)*ds(i) + (R(i)/2+(R(i-1)-R(i))/6)*ds(i-1));
                    end; 
                    A(1,1) = A(1,1) - dpdz*cos(alpha(1))*( (R(1)/2+(R(2)-R(1))/6)*ds(1) );
                    A(N,1) =        - dpdz*cos(alpha(1))*(  (R(N)/2+(R(N-1)-R(N))/6)*ds(N-1));
                end
                
               %A(N,N) = m.gamma*(-1/ds(N-1)*(R(N)+R(N-1))/2);
               KKi = m.gamma*(Ka(N)^2+Kb(N)^2)-dpdz*cos(alpha(N));  
               dmKKi =  m.gamma*(Ka(N-1)^2+Kb(N-1)^2)-dpdz*cos(alpha(N-1))-KKi;
               ri = R(N) ; dmri = R(N-1)-ri; 
               A(N,N) =  (20*KKi*ri+5*KKi*dmri+5*dmKKi*ri+2*dmKKi*dmri)*ds(1)/60;
               A(N,N-1) = A(N-1,N);
                
                
                switch(m.typeend)
                    
                    case('pined')  % "pinned" case
                    A(N,N) = -1;
                    A(N,N-1) = 0;
                    
                     case('axisX') %for sessile drop  : point 1 is on symmetry axis
                         % nothing to add
                    % A(N,N) = -1/ds(N-1)*(R(N-1)/2);  %%% facteur 4 car 1/R d eta / ds = d^2 eta / d s^2 (mais en fait 2...) 
                    % A(N,N-1) = A(N-1,N);
                    % A(N,N) = A(N,N)+ (   (m.gamma*(Ka(N)^2+Kb(N)^2)-dpdz*cos(alpha(N)))*R(N) ...
                    %                 + (m.gamma*(Ka(N-1)^2+Kb(N-1)^2)-dpdz*cos(alpha(N-1)))*R(N-1) )/2 * (ds(N-1)/3); 
                    
                    case('angle') % for hanging drop when imposing the angle at last point
        
                    %A(N,N-1) = m.gamma*(R(N)+R(N-1))/2*1/ds(N-1);
                    %A(N,N) = - m.gamma*R(N)*(1/ds(N-1));% 
                    %A(N,N-1) = A(N,N-1)  + (  (m.gamma*(Ka(N)^2+Kb(N)^2)-dpdz*cos(alpha(N)))*R(N) ...
                    %                        + (m.gamma*(Ka(N-1)^2+Kb(N-1)^2)-dpdz*cos(alpha(N-1)))*R(N-1) )/2 * (ds(N-1)/6);  
                    %A(N,N) = A(N,N) + (m.gamma*(Ka(N)^2+Kb(N)^2)- dpdz*cos(alpha(N)))*R(N)*(ds(N-1))/3;
                    %A(N,N-2) = 0;%pour etre sur
                    % terme supplementaire venant du terme de bord de l'IPP
                    % est a mettre dans le RHS
                    
                    % terme supplementaire traduisant la translation
                    % verticale pour recoller ?? la surface
                    %for i=2:N-1
                    %A(i,N) = A(i,N) - dpdz*cos(alpha(N))*R(i)*(ds(i-1)+ds(i))/2;
                    %end; 
                    %A(1,N) = A(1,N) - dpdz*cos(alpha(N))*(R(1)+R(2))/2*(ds(1))/2;
                    %A(N,N) = A(N,N) - dpdz*cos(alpha(N))*(R(N)+R(N-1))/2*(ds(N-1))/2;
                    
                     for i=2:N-1
                    A(i,N) = A(i,N) - dpdz*cos(alpha(N))*( (R(i)/2+(R(i+1)-R(i))/6)*ds(i) + (R(i)/2+(R(i-1)-R(i))/6)*ds(i-1));
                    end; 
                    A(1,N) = A(1,N) - dpdz*cos(alpha(N))*( (R(1)/2+(R(2)-R(1))/6)*ds(1) );
                    A(N,N) = A(N,N) - dpdz*cos(alpha(N))*(  (R(N)/2+(R(N-1)-R(N))/6)*ds(N-1));
                    
                end
                
                
                
                m.AP = A;
                m.AV(1:N,1:N) = A;
                
                % derniere ligne : 1/(2*pi) dV 
                for i=2:N-1
                    % m.AV(N+1,i) =  R(i)*(ds(i)+ds(i-1))/2; % test sign
                    % ca semble + logique d'avoir
                    % (R(i)+R(i-1)/2)*ds(i-1)/2 + (R(i)+R(i+1))/2*ds(i)/2  
                    % a tester... 
                    m.AV(N+1,i) = ds(i)*(R(i)/2+(R(i+1)-R(i))/6)+ds(i-1)*(R(i)/2+(R(i-1)-R(i))/6);
                end;
                m.AV(N+1,1) = ds(1)*(R(1)/2+(R(2)-R(1))/6); 
                m.AV(N+1,N) = ds(N-1)*(R(N)/2+(R(N-1)-R(N))/6); 
                
                %switch(m.typestart)
                %    case('axisX')
                %    m.AV(N+1,1) = ds(1)*(R(1)/2+(R(2)-R(1))/6); 
                %end
                
                %if(strcmpi(m.typeend,'pined')==1)
                %       m.AV(N+1,N) = R(N-1)/2*ds(N-1)/2; 
                %elseif (strcmpi(m.typeend,'axisX')==1)
                %    m.AV(N+1,N) = 0*R(N-1)/2*ds(N-1)/2; % test sign
                %else
                %    m.AV(N+1,N) = 0*R(N-1)/2*ds(N-1)/2;
                %end
                %m.AV(N+1,N+1) = 0;

                % last column : effet of dP
                for i=1:N
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
                    m.RHS(i) = m.gamma*(m.Ka(i)+m.Kb(i))-m.P+m.rhog*m.Z(i);
                end
                
                switch(m.typestart)
                    case('axisX')
                    m.RHS(1) = m.gamma*(m.Ka(1)+m.Kb(1))-m.P+m.rhog*m.Z(1);
                    case('pined')
                    m.RHS(1) = 0;
                    case('angle')
                    m.RHS(1)   =  (m.alpha(1)-m.beta); % a verifier
                end

                switch(m.typeend)
                    case('axisX')
                    m.RHS(m.N) = m.gamma*(m.Ka(m.N)+m.Kb(m.N))-m.P+m.rhog*m.Z(m.N);
                    case('pined')
                    m.RHS(m.N) = 0;
                    case('angle')
                    m.RHS(m.N)   =  (m.alpha(m.N)-m.beta);
                end

                % For V continuation : volume constraint
                if (m.typecont=='V')
                    m.RHS(m.N+1) = -(m.Vtarget-m.V);
                end   
                
                % For S continuation : arclength step
                if (m.typecont=='S')
                    m.RHS(m.N+1) = 0;
                end

                
                
                case('FE') %FE
                    
                for i=2:m.N-1
                YLi = (m.gamma*(m.Ka(i)+m.Kb(i))-m.P+dpdz*m.Z(i));
                dpYLi = (m.gamma*(m.Ka(i+1)+m.Kb(i+1))-m.P+dpdz*m.Z(i+1)) - YLi;
                dmYLi = (m.gamma*(m.Ka(i-1)+m.Kb(i-1))-m.P+dpdz*m.Z(i-1)) - YLi;
                ri = m.R(i) ; dpri = m.R(i+1)-ri ; dmri = m.R(i-1)-ri ;
                m.RHS(i) = (6*YLi*ri+2*YLi*dpri+2*dpYLi*ri+dpYLi*dpri)*m.ds(i)/60 ...
                                    + (6*YLi*ri+2*YLi*dmri+2*dmYLi*ri+dmYLi*dmri)*m.ds(i-1)/60 ;
                %    m.RHS(i) = (m.gamma*(m.Ka(i)+m.Kb(i))-m.P+dpdz*m.Z(i))*m.R(i)*(m.ds(i-1)+m.ds(i))/2; 
                end
                YLi = (m.gamma*(m.Ka(1)+m.Kb(1))-m.P+dpdz*m.Z(1));
                dpYLi = (m.gamma*(m.Ka(2)+m.Kb(2))-m.P+dpdz*m.Z(2)) - YLi;
                ri = m.R(1) ; dpri = m.R(2)-ri ; 
                m.RHS(1) = (6*YLi*ri+2*YLi*dpri+2*dpYLi*ri+dpYLi*dpri)*m.ds(1)/60;
                    
                YLi = (m.gamma*(m.Ka(m.N)+m.Kb(m.N))-m.P+dpdz*m.Z(m.N));
                dmYLi = (m.gamma*(m.Ka(m.N-1)+m.Kb(m.N-1))-m.P+dpdz*m.Z(m.N-1)) - YLi;
                ri = m.R(m.N) ; dmri = m.R(m.N-1)-ri ; 
                m.RHS(m.N) = (6*YLi*ri+2*YLi*dmri+2*dmYLi*ri+dmYLi*dmri)*m.ds(m.N-1)/60;
                
                switch(m.typestart)
                    case('axisX')
%                    m.RHS(1) = (m.gamma*(m.Ka(1)+m.Kb(1))-m.P+dpdz*m.Z(1))*((m.R(1)+m.R(2))/2)*(m.ds(1)/2);
                    case('pined')
                    m.RHS(1) = 0;
                    case('angle')
                    % a verifier
                    m.RHS(1) =m.RHS(1)+m.gamma*(m.alpha(1)-m.beta)*(m.R(1));
                end
                
                
                
                switch(m.typeend)
                    case('pined')
                    m.RHS(m.N) = 0;
                    case('angle')
                    % a verifier
                    m.RHS(m.N) = m.RHS(m.N) -m.gamma*(m.alpha(m.N)-m.beta)*(m.R(m.N));
                    case('axisX')
%                    m.RHS(m.N) =0;
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
                if(ismember(30,m.whichfigures)||m.idebug>49)
                figure(30);
                hold on;
                plot(m.P,real(m.lambdaminP),'r+',m.P,real(m.lambdaminV),'b+',m.P,real(m.lambdamin1),'g+');
                plot(m.P,imag(m.lambdaminP),'rx',m.P,imag(m.lambdaminV),'bx',m.P,imag(m.lambdamin1),'gx');
                end
            end
            
            % detection bifurcation P
            if (real(m.lambdaminP)*real(lambdaminPans)<0&&abs(m.lambdaminP-lambdaminPans)<1&&abs(imag(m.lambdaminP))<0.5) 
                VbV = (m.V*lambdaminPans - m.Vans*m.lambdaminP)/(lambdaminPans-m.lambdaminP);
                PbV = (m.P*lambdaminPans - m.Pans*m.lambdaminP)/(lambdaminPans-m.lambdaminP);
                if(ismember(20,m.whichfigures)||m.idebug>49)
                figure(20);
                hold on;
                plot(PbV,VbV,'xr','MarkerSize',10);
                end
                if(ismember(26,m.whichfigures)||m.idebug>49)
                figure(26);
                hold on;
                plot(PbV,VbV/m.Vref,'xr','MarkerSize',10);
                end
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
                
                if(ismember(20,m.whichfigures)||m.idebug>49)
                figure(20);
                hold on;
                plot(real(PbV),real(VbV),'+b','MarkerSize',10);
                end
                if(ismember(25,m.whichfigures)||m.idebug>49)
                figure(25);
                hold on;
                plot(real(PbV),real(VbV/m.Vref),'+b','MarkerSize',10);
                end
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
                if(ismember(20,m.whichfigures)||m.idebug>49)
                figure(20);
                hold on;
                plot(PbV,VbV,'*g','MarkerSize',10);
                end
                if(ismember(25,m.whichfigures)||m.idebug>49)
                figure(25);
                hold on;
                plot(PbV,VbV/m.Vref,'*g','MarkerSize',10);
                end
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
            if(ismember(numfig,m.whichfigures)||m.idebug>49)
            figure(numfig);
            plot(m.R,m.Z,plotstyle);
            hold on;
            grid on
            end
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
        
          function obliqueG(m,eps);

            % function to compute and plot the shape of the meniscus in a 
            figure(40);
            plot(m.R,m.Z,'k',-m.R,m.Z,'k');
            hold on;
            m.RHS = m.R*m.rhog;
            m.RHS(m.N)=0;
            m.eta(1:m.N) = (m.A1(1:m.N,1:m.N)\m.RHS(1:m.N)')';
            plot( m.R+eps*m.eta.*sin(m.alpha), m.Z-eps*m.eta.*cos(m.alpha),'r');
            plot(-m.R+eps*m.eta.*sin(m.alpha), m.Z+eps*m.eta.*cos(m.alpha),'r');

        end %function obliqueG
         
        function resetfigs(m);
            % to reset figures
 
            close all;
            
         if(ismember(10,m.whichfigures)||m.idebug>49)
            figure(10);
            title('Menisci shapes');
            xlabel('R');
            ylabel('Z');
%            axis equal;
            hold on;
         end
            
         if(ismember(11,m.whichfigures)||m.idebug>49)
            figure(11);
            title('Menisci shapes (all computed)');
            xlabel('R');
            ylabel('Z');
            hold on;
         end
             
         if(ismember(20,m.whichfigures)||m.idebug>49)
            figure(20);
            title('P/V curves');
            xlabel('P');
            ylabel('V');
            hold on;
         end
    
        if(ismember(21,m.whichfigures)||m.idebug>49)
            figure(21);
            title('P/V results (dots)');
            xlabel('P');
            ylabel('V');
            hold on;
        end
        
        if(ismember(22,m.whichfigures)||m.idebug>49)
            figure(22);
            title('P/H results');
            xlabel('P');
            ylabel('H');
            hold on;
        end
        
        if(ismember(23,m.whichfigures)||m.idebug>49)   
            figure(23);
            title('P/R0 results');
            xlabel('P');
            ylabel('R0');
            hold on;
        end
        
        if(ismember(25,m.whichfigures)||m.idebug>49)
            figure(25);
            title('bif. points in PV diagram');
            xlabel('P');
            ylabel('V');
            hold on;
        end
        
        if(ismember(26,m.whichfigures)||m.idebug>49)  
            figure(26);
            title('P/V* results');
            xlabel('P');
            ylabel('V*');
            hold on;
        end
        
         if(ismember(28,m.whichfigures)||m.idebug>49)   
            figure(28);
            title('H/R0 results');
            xlabel('R0 (radius of droplet)');
            ylabel('H (height of droplet)');
            hold on;
        end
        
        if(ismember(30,m.whichfigures)||m.idebug>49)
            figure(30);
            title('Eigenvalues of matrices AP, AV, A1');
            ylabel('lambdaP, lambdaV, lambda1');
            xlabel('P');
            hold on;
        end
       
         if(ismember(121,m.whichfigures)||m.idebug>49)
            figure(121);
            title('Static angle vs volume (for droplet with pinned contact line)');
            ylabel('\theta_s (deg)');
            xlabel('V');
            hold on;
         end
       
  
            
        end % function resetfigs
           
    end % methods

end % class meniscus 