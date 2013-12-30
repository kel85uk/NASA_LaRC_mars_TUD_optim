function ydot=sysSolveCRTBP(t,y,options,flag, mu)

%Not needed but lets you recycle more code...
G=1;

%This file contains the right hand side for the system which defines the
%state transition matrix for the spatial CRTBP.  

%-------------------------------------------------------------------------
%The structure of the system is 
%
%          A'=Df*A, x'=f(x)    
%
%Where A,Df, are 6X6 matrices and A' is the derivative of A.
%Df is the derivative of the N-body vector field f and has the form
%
%                    |  0   I  |  
%               Df = |         |
%                    | G    K  | 
%
%where all submatrices are 3X3 and G has depends only on the position
%vectors of the bodies (all be it in a complicated way), and x'=f(x)
% is the spatial CRTBP.  All of this is combined onto one 
%system of 1st order ODEs. The right hand side for the system is coded 
%in the remainder of the file
%------------------------------------------------------------------------


%FIRST
%The matrix G is computed
%This matrix depends on the positions of the N-body problem.
%These positions are contained in the 6 entries y(145:150).
    
%put the positions aside
x(1:3)=y(37:39);
   
%Now compute the matrix G.  Since 'G' already denotes the gravatational 
%constant call the matrix G 'GMatrix'. 
%This is done by calling 'G_CRTBP.m'
GMatrix=G_CRTBP(x, mu);

%SECOND
%The right hand side for the state transition system is computed
%To do this construct a matrices 'A' 'I' 'K' and 'O', where 'A' contains the 
%variables and 'O' is a 3X3 matrix of zeros and 'I' is the 3X3 identity.   
O=zeros(3);
I=eye(3);
K=[0, 2, 0;
  -2, 0, 0;
   0, 0, 0];


%Now the complete Jacobian of f is assembled
Df=0;
Df=[ O,       I;
     GMatrix, K];

%Make A
for i=1:6
    for j=1:6
        A(i,j)=y(6*(i-1)+j);
    end
end

%Then compute the 6X6 matrix Df*A, is named DfA.
DfA=0;
DfA=Df*A;

%This has to be put into vector format. We temporaly place the results in 
%the 36-vector 'a'.  Later this will be the first 36 components of ydot.
a=0;
for i=1:6
    for j=1:6
        a(6*(i-1)+j)=DfA(i,j);
    end
end    

%THIRD
%The last 6 entries are the vector field for the CRTBP,
%These are stored in 'c'.



%-------------------------------------------------------------------------
%--------The following code is essentialy the code from------------------
%-------'CRTBP.m' copied here and adjusted for 2 dimensions----------------
%-------------------------------------------------------------------------

%the distances
r1=sqrt((mu+y(37))^2+(y(38))^2+(y(39))^2);
r2=sqrt((1-mu-y(37))^2+(y(38))^2+(y(39))^2);
%masses
m1=1-mu;
m2=mu;


%constructs a vector whose first 3 entries are the velocities and whose
%last 3 entries are the accelerations

c=[y(40); 
    y(41); 
    y(42); 
    y(37)+2*y(41)+G*m1*(-mu-y(37))/(r1^3)+G*m2*(1-mu-y(37))/(r2^3); 
    y(38)-2*y(40)-G*m1*(y(38))/(r1^3)-G*m2*y(38)/(r2^3); 
    -G*m1*y(39)/(r1^3)-G*m2*y(39)/(r2^3)];

%------------------------------------------------------------------------
%-----------------------------'CRTBP.m'----------------------------------
%------------------------------------------------------------------------

%Put it all toghether and pass back to integrator

ydot=[a';c];