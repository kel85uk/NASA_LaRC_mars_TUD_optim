function[a,dVtot,dV0,dV3] = Lambert_targ(dtheta,deltat,mu1,mu2,r1,r2,a1_park,a2_park,e1_park,e2_park)
% function[a,dVtot,dV0,dV3] = Lambert_targ(dtheta,deltat,mu1,mu2,r1,r2,a1_park,a2_park,e1_park,e2_park)
%
% This function will compute a solution to an arbitrary lambert targeting
% problem, assuming coplanar departure and arrival locations, and a
% coplanar transfer. It is also assumed that the 
%
% Inputs are as follows:
%
% r1,r2 = orbital radius departure/arrival planet [AU]
% dtheta = angular separation departure and arrival points [deg]
% deltat = time of flight in days [solar days = 86400 s]
% a1_park,a2_park = semi-major axis of parking orbit at departure/arrival planet [km]
% e1_park,e2_park = eccentricity of parking orbit at departure/arrival planet
% mu1,mu2 = gravitational parameter of departure/arrival planet [km3/s2]
%
% Written by Harry Linskens
%
% VERSION HISTORY
% V1.0: Written 15-12-2012
% V1.1: Updated to include deltaV calculations, 18-12-2012.
% V1.2: Modified to accept elliptical parking orbits, 16-05-2013
%--------------------------------------------------------------------------
format longg

% Initialize parameters
AU = 149.597871e6;                  % [km]
mu_sun = 1.3271244e11;              % [km3/s2]
r1 = r1*AU;                         % [km]
r2 = r2*AU;                         % [km]
dtheta = deg2rad(dtheta);           % [rad]
daylength = 86400;                  % [s]
deltat = deltat*daylength;          % [s]
x0 = 0;

% The chord between the departure and arrival points must now be
% determined. This is done using the cosine rule. Then, deltat is
% normalized.
c = sqrt( r1^2 + r2^2 - 2*r1*r2*cos(dtheta) );      % [km]
s = (r1 + r2 + c)/2;                                % [km]
T = sqrt(8*mu_sun/(s^3))*deltat;



% Further parameters can then be defined. Also, some values for the loop
% are initialized.
q = (sqrt(r1*r2)/s) * cos(dtheta/2);
F = 1;
tolerance = 1e-3;
diffstep = 1e-3;
counter = 1;

%--------------------------------------------------------------------------
% Initiate a first guess of x0
x = x0;
    E_lam = x^2-1;
    y = sqrt(abs(E_lam));
    z = sqrt(1-q^2 + q^2*x^2);
    f = y*(z - q*x);
    g = x*z - q*E_lam;

    % Determine value of d
    if E_lam < 0
        d = atan2(f,g);
    elseif E_lam > 0
        d = real(log(f+g));
%     else
%         error('lmbrt_targ_E','E_lam is equal to zero, and the notes dont cover this')
    end
T0 = 2*(x-q*z-(d/y))/E_lam;

if T0 > T
    x0 = T0*(T0-T)/(4*T);
elseif T0 < T
    x0 = -(T-T0)/(T-T0+4);
elseif T0 == T
    x0 = 0;
end

%--------------------------------------------------------------------------
% Iterate using Halley's method
x = x0;
while abs(F) > 0+tolerance
    E_lam = x^2-1;
    y = sqrt(abs(E_lam));
    z = sqrt(1-q^2 + q^2*x^2);
    f = y*(z - q*x);
    g = x*z - q*E_lam;
    
    % V1.3: Hard-coded this in order to prevent singularities when
    % assessing near-parabolas.
    if z > 1-0.0001
        tolerance = 0.005;
    end

    % Determine value of d
    if E_lam < 0
        d = atan2(f,g);
    elseif E_lam > 0
        d = real(log(f+g));
%     else
%         error('lmbrt_targ_E','E_lam is equal to zero, and the notes dont cover this')
    end
    % F can then be defined
    F = T - ((2*(x - q*z - (d/y))) / (E_lam));
    
    % Halleys method requires the first and second derivatives of F. This
    % is done numerically. First, the first derivative
    xprime = x+diffstep;
    E_lamprime = xprime^2-1;
    yprime = sqrt(abs(E_lamprime));
    zprime = sqrt(1-q^2 + q^2*xprime^2);
    fprime = yprime*(zprime - q*xprime);
    gprime = xprime*zprime - q*E_lamprime;
    % Determine value of dprime
    if E_lamprime < 0
        dprime = atan2(fprime,gprime);
    elseif E_lam > 0
        dprime = real(log(fprime+gprime));
%     else
%         error('E_lam is equal to zero, and the notes dont cover this')
    end
    Fprime = T - ((2*(xprime - q*zprime - (dprime/yprime))) / (E_lamprime));
    dFdx = (Fprime-F) / diffstep;
    
    % Now, the second derivative
    xdprime = xprime+diffstep;
    E_lamdprime = xdprime^2-1;
    ydprime = sqrt(abs(E_lamdprime));
    zdprime = sqrt(1-q^2 + q^2*xdprime^2);
    fdprime = ydprime*(zdprime - q*xdprime);
    gdprime = xdprime*zdprime - q*E_lamdprime;
    % Determine value of ddprime
    if E_lamdprime < 0
        ddprime = atan2(fdprime,gdprime);
    elseif E_lam > 0
        ddprime = real(log(fdprime+gdprime));
%     else
%         error('E_lam is equal to zero, and the notes dont cover this')
    end
    Fdprime = T - ((2*(xdprime - q*zdprime - (ddprime/ydprime))) / (E_lamdprime));
    d2Fdx2 = (Fdprime - Fprime) / diffstep;
    
    % Finally, Halleys method can be applied to determine the next value of
    % x
    xnew = x - (F*dFdx/(2*dFdx^2 - F*d2Fdx2));
%     xnew = x - F/dFdx;
    x = xnew;
    counter = counter + 1;
end

x_final = x;

% With x_final found, the orbital parameters and inertial velocities can be
% determined.
a = ((s/2)/(1-x_final^2))/AU;        % [AU]

% Interim parameters are defined.
gamma = sqrt(mu_sun*s/2);
rho = (r1-r2)/c;
sigma = 2*sqrt(r1*r2/c^2)*sin(dtheta/2);

Vr1 = gamma*( (q*z-x) - rho*(q*z + x)) / r1;
Vr2 = -gamma*( (q*z-x) + rho*(q*z + x)) / r2;
Vt1 = gamma*sigma*(z+q*x)/r1;
Vt2 = gamma*sigma*(z+q*x)/r2;


%% HELIOCENTRIC SCALE
% The velocities of the planets in orbit can be determined.
Vrp = 0;
Vtp1 = sqrt(mu_sun/r1);
Vtp2 = sqrt(mu_sun/r2);

% The hyperbolic excess velocities can then be computed.
Vinf1 = sqrt((Vr1-Vrp)^2 + (Vt1-Vtp1)^2);
Vinf2 = sqrt((Vr2-Vrp)^2 + (Vt2-Vtp2)^2);


%% PLANETOCENTRIC SCALE - DEPARTURE
% Assume the necessary burns occur at perigee. Then, compute escape
% velocity and hyperbolic perigee velocity.
Vesc1 = sqrt(2*mu1/(a1_park*(1-e1_park)));                      % [km/s]
Vper1 = sqrt(mu1* (2/(a1_park*(1-e1_park)) - 1/a1_park) );      % [km/s]
Vhyp1 = sqrt(Vesc1^2 + Vinf1^2);                                % [km/s]

% Now, the first maneuver delta-V can be determined
dV0 = Vhyp1 - Vper1;                                   % [km/s]


%% PLANETOCENTRIC SCALE - ARRIVAL
% Assume the necessary burns occur at perigee. Then, compute escape
% velocity and hyperbolic perigee velocity.
Vesc2 = sqrt(2*mu2/(a2_park*(1-e2_park)));                      % [km/s]
Vper2 = sqrt(mu2* (2/(a2_park*(1-e2_park)) - 1/a2_park) );      % [km/s]
Vhyp2 = sqrt(Vesc2^2 + Vinf2^2);                                % [km/s]

% Now, the final maneuver delta-V can be determined
dV3 = Vhyp2 - Vper2;                                   % [km/s]


%% FINAL PROCESSING
dVtot = dV0 + dV3;                           % [km/s]


