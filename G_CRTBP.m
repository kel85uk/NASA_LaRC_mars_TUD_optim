function GMatrix=G_CRTBP(x, mu)

%function returns the matrix 'G' for the CRTBP  
%with parameter mu.


%the distances
r1=sqrt((x(1)+mu)^2+x(2)^2+x(3)^2);
r2=sqrt((x(1)-(1-mu))^2+x(2)^2+x(3)^2);

%If U is the potential for the CRTBP then this code is computing 
%the derivative of the gradient of U.  The gradient has three 
%components which we will call u1, u2, u3.  The differential of the 
%gradient is the matrix of partials of these functions.  These will 
%be denoted by u1_x, u1_y, and so forth.  

u1_x=1-(1-mu)*(1/(r1^3)-3*((x(1)+mu)^2)/(r1^5))-mu*(1/(r2^3)-3*((x(1)-(1-mu))^2)/(r2^5));

u2_y=1-(1-mu)*(1/(r1)^3-3*x(2)^2/r1^5)-mu*(1/r2^3-3*x(2)^2/r2^5);

u3_z=(-1)*(1-mu)*(1/(r1)^3-3*x(3)^2/r1^5)-mu*(1/r2^3-3*x(3)^2/r2^5);

u1_y=3*(1-mu)*x(2)*(x(1)+mu)/r1^5+3*mu*x(2)*(x(1)-(1-mu))/r2^5;

u1_z=3*(1-mu)*x(3)*(x(1)+mu)/r1^5+3*mu*x(3)*(x(1)-(1-mu))/r2^5;

u2_z=3*(1-mu)*x(2)*x(3)/r1^5+3*mu*x(2)*x(3)/r2^5;


%equality of mixed partials gives (as all the terms are already partials
%of the potential function);
u3_y=u2_z;

u2_x=u1_y;

u3_x=u1_z;



%Then (as mentioned) G is the matrix of partials

GMatrix=[u1_x, u1_y, u1_z;
         u2_x, u2_y, u2_z;
         u3_x, u3_y, u3_z];
     
         
         
         
         
         