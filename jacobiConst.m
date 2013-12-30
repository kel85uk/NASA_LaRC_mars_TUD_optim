function C=jacobiConst(y,v,mu)

%File computes Jacobi Energy for a given state (x,v)

%the distances
r1=sqrt((mu+y(1,1))^2+(y(2,1))^2+(y(3,1))^2);
r2=sqrt((y(1,1)-(1-mu))^2+(y(2,1))^2+(y(3,1))^2);

%Compute the Jacobi Energy
C=-(v(1,1)^2 + v(2,1)^2+v(3,1)^2)/2 + 2*((y(1,1)^2 + y(2,1)^2)/2 + (1-mu)/r1 + mu/r2);

