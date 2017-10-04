
classdef meniscus
    
    properties
        Ndim = 3
        discretization = 'FD'
        dpdz = 1
        typestart = 'axis' 
        typeend = 'pinned'
        N
        R
        Z
        P 
        V
        alpha 
        Ka
        Kb
        s0
        ds
        AV
        AP
        A1
        betastart
        betaend
        paramscan = 'Bo'
        conv = 1;
        
    end
    
    
    methods
        function m = meniscus(N,Rdisk) % constructor for a plane, circular interface
            m.N = N;
            t = 0:1/N:1;
            m.R = Rdisk*t;
            m.Z = 0*t;
            m.V = 0;
            m.P = 0;
        end
           
        function m = calcgeom(m);
           ds(1) = sqrt((R(2)-R(1))^2+(Z(2)-Z(1))^2);
            s0(1) = 0;
            
            for i=2:N-1
                ds(i) = sqrt( (R(i+1)-R(i))^2+(Z(i+1)-Z(i))^2);
                s0(i) = s0(i-1)+ds(i-1);
            end

            s0(N) = s0(N-1)+ds(N-1);
     

    
            for i=1:N-1
                alphaS(i) =  atan2(Z(i+1)-Z(i),R(i+1)-R(i)); % angle at the middle of the side

            end       
 
            unwrap(alphaS); 

    for i=2:N-1
        alpha(i) =  (ds(i)*alphaS(i-1) +  ds(i-1)*alphaS(i))/(ds(i-1)+ds(i));
        Ka(i) = (alphaS(i)-alphaS(i-1))/((ds(i)+ds(i-1))/2);
        if (nbdim == 3) 
            Kb(i) = sin(alpha(i))/R(i);
        else
            Kb(i) = 0;
        end
    end;
    
    % i = N : pinned point (geom. quantities only needed for visualization and for angle-fixed case)
    
  %  alpha(N) =2*alpha(N-1)-alpha(N-2)
    alpha(N) = alpha(N-2)+2*ds(N-1)*(Ka1(N-1)+Kb1(N-1)-sin(alpha(N-1))/R(N-1));
    Ka(N) =2*Ka1(N-1)-Ka1(N-2); 
    Kb(N) =2*Kb1(N-1)-Kb1(N-2); 
    
  %  alpha(N-2:N)
    
    if(abs(R(1))>1e-6)
    %  i = 1 : pinned point
    alpha(1) = 2*alpha(2)-alpha(3); 
    Ka(1) = 2*Ka1(2)-Ka1(3); 
    Kb(1) = 2*Kb1(2)-Kb1(3); 

    else
     %  i = 1 : axis   
     alpha(1) = 0.;
    Ka(1) = 2*atan2(Z(2)-Z(1),R(2)-R(1) )/ds(1); % facteur 2 rajouté 
      if (nbdim == 3) 
          Kb(1) = Ka1(1);
      end
    end     
    
    unwrap(alpha);
        
        end
    
        function m = Newton_P(mans,P);
            m = mans;
            m.P = P;
            m = calcgeom(m);
            m.conv = 1;
            if(m.conv == 0)
                m = mans;
                m.conv = 0;
                disp ('divegence !')
            end
        end
    end
end