function xi = calcxi(R,Z,xi0,xiN);
[alpha,ds,s0,Ka,Kb] = calcgeom(R,Z);

N = length(R);
curvelength = s0(N);

%figure(27);
%plot(R,Z,'bo')
%hold on;



for i=1:N
    xi(i)=xi0+(i-1)/(N-1)*(curvelength+xiN)-s0(i);
end


ds(N) = ds(N-1);

%figure(456);
%
%plot(s0,ds,'r',s0,xi,'b')
%figure(457);
%plot(s0);

%plot(R,Z,'r+')

