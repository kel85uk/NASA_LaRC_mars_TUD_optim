function A=stateTransCRTBP(t0, tf, x, mu)

%Compute the state transition matrix at time tf for the 
%path x(t) with x(t0)=x0

tspan=[t0,tf];                        %time span over which to run the integration

%--------------------------------------------------------------------------
%The state transition matrix is determined by a 6X6 matrix ODE
%This gives rise first to a system of 36 ODEs.  However the matrix
%ODE is nonautonomous and depends on a particular solution (trajectory)
%of the CRTBP.  This trajectory is itself the solution of a system
%of 6 ODEs.  The two systems are solved simultaneously, giving an
%autonomous system of 36+6=42 ODEs.
%--------------------------------------------------------------------------

%set up the initial condition, y0 for the 42 component system

y0=0;

%since the initial condition for the 6X6 matrix system is the 6X6
%identity matrix the initial condition vector is very sparse.  The first
%36 components are 0s and 1s, and the last 6 are the initial conditions
%from the CRTBP.

 %Initialize a 6X6 identity matrix
I=eye(6);                           

%put the entries, row by row, into y0.
for i=1:6
   for j=1:6;
      y(6*(i-1)+j)=I(i,j); 
   end    
end    

%the inital conditions for the particular orbit of the CRTBP have
%been passed in as 'x'. This to the end of y.

y(37:42)=x';

%Convert to a column vector and pass the initial conditions to the
%integrator
y0=y';                                                     %inital condition
options=odeset('RelTol',1e-13,'AbsTol',1e-22);             %set tolerences
[t,Y]=ode113('sysSolveCRTBP',tspan,y0,options,[], mu);   %integrate the system

%[t,Y] is a huge matrix.  It is made up of a row for each time and 42
%columns.  But the desired data is the state transition matrix at the final
%time.  Then the first 36 entries of the last row of Y must be put into a
%6X6 matrix which will be passed back to the caller

%find out the row number of tf
b=size(Y);        %returns the number of rows and columns as a row vector
m=b(1,1);         %so 'm' is the number of rows
c=Y(m,1:36);     %so 'c' is a vector containing the first 36 entries of the last row of Y 

%now c has to be made into a 12X12 matrix
d=0;
for i=1:6
    for j=1:6
        d(i,j)=c(6*(i-1)+j);
    end    
end

%d is the state transition matrix at the final time and is passed back to
%the caller
A=d;
    

