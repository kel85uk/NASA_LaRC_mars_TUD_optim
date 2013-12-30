function ydot=CRTBP(t,y,options,flag,G,mu)

%the distances
r1=sqrt((mu+y(1))^2+(y(2))^2+(y(3))^2);
r2=sqrt((1-mu-y(1))^2+(y(2))^2+(y(3))^2);
%masses
m1=1-mu;
m2=mu;

ydot=[y(4); 
    y(5); 
    y(6); 
    y(1)+2*y(5)+G*m1*(-mu-y(1))/(r1^3)+G*m2*(1-mu-y(1))/(r2^3); 
    y(2)-2*y(4)-G*m1*(y(2))/(r1^3)-G*m2*y(2)/(r2^3); 
    -G*m1*y(3)/(r1^3)-G*m2*y(3)/(r2^3)];