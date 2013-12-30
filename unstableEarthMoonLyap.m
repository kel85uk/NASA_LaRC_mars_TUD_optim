%given an initial position of a periodic orbit this program
%computes stable and stable manifolds of the orbit (assuming they exist...)

%clear the workspace
clear;
%output long
format long;

%--------------------------------------------------------------------------
%----------------------NOTES-----------------------------------------------
%--------------------------------------------------------------------------

%NAME:  unstableEarthMoonLyap.m

%DEPENDANCIES: The program calls 'CRTBP.m', 'librationPoints.m', 
%'stateTransCRTBP.m'and 'jacobiConst.m'.  The program 'stateTransCRTBP
%calls 'sysSolveCRTBP', which calls 'G_CRTBP.m'.

%These must be in the MatLab search path in order for the present
%program to run.

%USE: Once the program and it's dependencies are in the search path simply
%type:
%                  unstableEarthMoonLyap
%
%from the matlab The program also keeps all the
%data points for the orbit and it's unstable manifold.  This can be
%obtained by typing 'who' at the command prompt after the program has run.command line.

%The program can be run with the default data currently set, or the uset
%can specify the initial condition of an unstable periodic orbit and it's
%period as described below.


%OUTPUT: The output is a plot of the Lyapunov orbit, it's unstable manifold 
%in the vacinity of the moon (secondary body), and the zero velocity curve
%at the energy level of the Lyapunov orbit. 


%PURPOSE: to compute the unstable manifold of a Lyapunov orbit in the
%planar Circular Restricted Three Body Problem.  This is mainly a
%visulation program as the points making up the unstable manifold is 
%plottednot but not kept.  The program can be easily adapted to the study
%of any orbit which has a single stable, and single unstable direction, and
%hence two dimensional stable and unstable manifolds.

%NOTE:  The word "halo" appears in several variable names throught the
%program.  This is because the program was originally developed to compute
%the stable and unstable manifolds of halo orbits.  However it is not
%necessary that orbit under consideration be a halo orbit.  In fact the
%default setting of the program are for a Lyapunov orbit near L1.

%--------------------------------------------------------------------------
%----------------------USER SPECIFIED--------------------------------------
%------------------------PARAMETERS----------------------------------------
%--------------------------------------------------------------------------



    %----------------PRIMARY CONTROL VARIABLES-----------------------------
    %how far to move in the unstable direction away from the orbit before
    %integrating;
    epsilon=2*10^(-8);

    %NOTE: The next data really decides what you are computing.  The user
    %should input 'r0' and 'v0' which are the initial position and velocity
    %for the periodic orbit they are interested in. The user also inputs the 
    %period of the orbit.  The defaults are for a Lyapunov orbit about L1
    %in the earth/moon system (this problem is planar). 

    %The initial condition for the halo orbit
    r0=[ 0.83946302646687;
                        0;
                        0];

    v0=[            0;
    -0.02596831282986;
                    0];

    %The period of the halo orbit is know as well (this should not be
    %changed unless the orbit is changed.
    haloPeriod= 2.69239959528586;
    
    
    %Parameterize the Manifold:  These two variables are the manifold
    %coordinates.  How long along the unstable holonomy and from how many 
    %and what starting points along the orbit
    
    
    %how long to move along the unstable orbit (how local
    %is the computation of the manifold)
    T=3.1*haloPeriod;

    %'k' is how many points to use in order to parameterize the periodic 
    %orbit.
    k=40;
    %----------------------------------------------------------------------


%--------------------------------------------------------------------------    
%----------------------SECONDARY PARAMETERS--------------------------------
%--------------------------------------------------------------------------    
 

%for the earth moon system the value of mu is approx 0.012277471
mu=0.012277471;

%Compute the Jacobi constant for the given initial conditions
C=jacobiConst(r0,v0,mu)

%Gravational constant
G=1;

%We are given data for the Sun/Earth/Moon system where the sun is the
%primary, and the earth and moon are combined into one body at their center
%of mass.  The earth moon system is refered to as the earth.  Then

GM_sun=1.327*10^(11);    
GM_earth=4.053*10^(5);
au=1.496*10^8;

%NOTE:  The previous paramters are not used in the present version of the
%program, but they are available in case one wished to convert into metric
%units at some point...

%The definition of mu gives
%mu=GM_earth/(GM_sun+GM_earth)



m1=1-mu;
m2=mu;



%Once 'G' and 'mu' are defined we can compute the location of the libration
%points and their Jacobi Integrals.  These will be handy when we set the
%value of the Jacobi Integral for the simulation below.

    %AXULARY CONSTANTS:
    %Compute the location of the five libration points as points in three
    %space.  These can be used to set the value of the Jacobi Constant.
    [L1, L2, L3, L4, L5]=librationPoints(mu);

    %Output location to the screen (if you want);
    outL1=L1;
    outL2=L2;
    outL3=L3;
    outL4=L4;
    outL5=L5;
    %Now compute the Jacobi Energy at each collinear libration point

        %The velocity at any equilibrium is zero.
        v_equil=[0; 0; 0];

    %Compute the energy at each libration point
    CL3=jacobiConst(L3, v_equil, mu);
    CL1=jacobiConst(L1, v_equil, mu);
    CL2=jacobiConst(L2, v_equil, mu);
    CL4=jacobiConst(L4, v_equil, mu);
    CL5=CL4 %so no extra computation is needed

%--------------------------------------------------------------------------
%-------------------COMPUTATIONS-------------------------------------------
%--------------------------------------------------------------------------



%Compute the periodic orbit
t0=0;
numSteps=300;
tspan=linspace(0,haloPeriod+0.1*haloPeriod,numSteps);
%numerical integration
y0=[r0; v0]';                                      %inital condition
options=odeset('RelTol',1e-13,'AbsTol',1e-22);     %set tolerences
[t,x_halo]=ode113('CRTBP',tspan,y0,options,[],G,mu);

%compute the monodromy matrix
Phi=stateTransCRTBP(t0, haloPeriod, [r0; v0], mu);

%compute the eigenvalues and eigenvectors

    %Get eigenvalues and eigenvectors
    [Phi_eigVecs, Phi_eigValOnDiag]=eigs(Phi);
    PhiEigs=diag(Phi_eigValOnDiag);


%NOTE:
%the first eigen value corrosponds to the unstable manifold
%the sixth (last) corrosponds to the stable direction. This is because
%MatLab orginizes the eigenvalues and their eigen vectore by the magnitude
%of the eigen values.  Then the first eigen value is the largest and the
%last is the smallest.  This assumes that the orbit has two dimensional 
%stable and unstable manifoldsl.  


figure
hold on

%get the parameterization points of the periodic orbit
h=round((numSteps-1)/k)
for n=1:k
    haloPoint(n,1:6)=x_halo((n-1)*h+1,1:6);
    haloTimes(n)=t((n-1)*h+1);
end




    I=eye(6);
    %The un-stable eigenvector is
    deltaU=Phi_eigVecs(1:6,1);

    %Compute an orbit on the unstable manifold
    %begining at the origin of the orbit
    xU_p(1,1:6)=haloPoint(1,1:6)'+epsilon*I*deltaU;


    %compute the unstable orbit
    %the positive ones;
    t0=0;
    numSteps=2000;
    tspan=linspace(0, T,numSteps);
    %numerical integration
    y0=xU_p(1,:)';                                      %inital condition
    options=odeset('RelTol',2.5e-13,'AbsTol',1e-22);     %set tolerences
    [t,xU_p1]=ode113('CRTBP',tspan,y0,options,[],G,mu);




    for n=2:k
        %compute the initial conditions for the manifold orbits
        initialconditions=n
        perturbationVector=stateTransCRTBP(t0, haloTimes(n), [r0; v0], mu)*deltaU;
        u_pV=perturbationVector/norm(perturbationVector);
        xU_p(n,1:6)=haloPoint(n,1:6)'+epsilon*u_pV;
    end



    %Compute and plot the unstable manifold

    for n=2:k
        manifoldLoop=n
        %the positive ones;
        t0=0;
        numSteps=2000;
        tspan=linspace(0, T,numSteps);
        %numerical integration
        y0=xU_p(n,1:6)';                                      %inital condition
        options=odeset('RelTol',2.5e-13,'AbsTol',1e-22);     %set tolerences
        [t,xU_pOrb]=ode113('CRTBP',tspan,y0,options,[],G,mu);

        %plot the unstable pos orbits as red lines
        plot3(xU_pOrb(:,1), xU_pOrb(:,2), xU_pOrb(:,3), 'r')
    end


%--------------------------------------------------------------------------
%---------------------Plotting---------------------------------------------
%--------------------------------------------------------------------------


%Contour Plot of the zero velocity surface
   [X,Y]=meshgrid(-1.1:0.002:1.2);
    r1=((X + mu).^2+Y.^2).^(1/2);
    r2=((X + mu-1).^2+Y.^2).^(1/2);
    Z=2*[(1/2)*(X.^2 + Y.^2)+(1-mu)./r1+mu./r2];
    contour(X,Y,Z,[C, C],'r')




%plot the libration points
plot3(L3(1), L3(2), L3(3), 'k*', L1(1),  L1(2), L1(3), 'k*')
plot3(L2(1), L2(2), L2(3), 'k*', L4(1), L4(2), L4(3), 'k*')
plot3( L5(1), L5(2), L5(3), 'k*')

%plot the planets/primaries 
plot3(1-mu, 0, 0,'g*')
plot3(-mu, 0, 0, 'b*')

%plot the periodic orbit in black dots
plot3(x_halo(:,1), x_halo(:,2), x_halo(:,3), 'k.')


