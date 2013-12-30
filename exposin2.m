function[r,gamma,TOF_mat,theta,dV_mat,a_T] = exposin2(r1,r2,Psi,TOF_req,timestep,N,k2,mu,n,gamma1)
% function[output] = exposin(r1,r2,t,N)
%
% This function will perform a Lambert-like targeting for LTP applications
% using exposins.
%
% r1,r2 =       orbital radius of departure,arrival planet [km]
% Psi =         angular separation of planets [deg]
% TOF_req =     required time of flight [sidereal days]
% timestep =    integration time step [sidereal days]
% N =           number of complete revolutions in transfer
% k2 =          assumed value of k2
% evaluations = number of values for gamma1 (excluding gamma1_min and
%               gamma1_max, so 1 for this will give one intermediate value
%               between the two limits). REMOVED IN V1.2
% mu =          gravitational parameter main body [km3/s2]
% n =           -1 for deceleration and 1 for acceleration
% gamma1 =      NOT column vector containing the initial flight path angles
%               [rad]
%
% [r,~,TOF_mat,theta,dV_mat,a_T] = exposin2(150000000,225000000,180,5*365,5,2,0.05,1.327e11,1,[0:0.05:0.0]')
%
% TEEEEEST! [r,gamma,TOF_mat,theta,dV_mat,a_T] = exposin2(1,2,1,1000,0.5,1,2/9,1.327e11,1,1
%
%
% Written by Harry Linskens
%
% VERSION HISTORY:
% V1.0: Written 24-12-2012. Spent 4 hours getting the geometry aspect to
%       work with a predefined theta-vector.
% V1.1: Written 26-12-2012. Added dynamical definition of theta using
%       thetadot. Also added TOF limit in integration loop.
% V1.2: Written 15-02-2013. Added acceleration and delta-V determination.
%       Saved as exposin2
%--------------------------------------------------------------------------
format longg

if n ~= 1 && n ~= -1
    error('n must be 1 or -1')
end

% With assumptions for value of k2 and gamma, some variables can be
% defined. Others can be initialized.
theta0 = 0; 
Psi = deg2rad(Psi);                     % [rad]

% TEEEEEST!
thetabar = Psi + 2*pi*N;                % [rad]
% thetabar = deg2rad(800);

daylength = 86164.1004;                 % [s]. sidereal day


% TEEEEEST!
% The constraint equation can be applied to limit options for gamma1.
Delta = 2*(1-cos(k2*thetabar))/(k2^4) - (log(r1/r2)^2);
gamma1_mintot = atan( (k2/2) * (-log(r1/r2)*cot(k2*thetabar/2) - sqrt(Delta)) );
gamma1_maxtot = atan( (k2/2) * (-log(r1/r2)*cot(k2*thetabar/2) + sqrt(Delta)) );

gamma1_min = gamma1(1,1);
gamma1_max = gamma1(end,1);

% EDIT V1.2, changed gamma1 to be an input
if gamma1_min < gamma1_mintot
    error('Minimum flight path angle is too small. Please increase.')
elseif gamma1_max > gamma1_maxtot
    error('Maximum flight path angle is too large. Please decrease.')
end


% TEEEEEST!
% Determine k1 and sign thereof.
k1 = sqrt( ((log(r1/r2)+tan(gamma1)*sin(k2*thetabar)/k2)/(1-cos(k2*thetabar))).^2 + (tan(gamma1).^2)/(k2^2) );
signk1 = (log(r1/r2)+tan(gamma1)*sin(k2*thetabar)/k2);

for i = 1:length(k1)
    if signk1(i,1) < 0
        k1(i,1) = -k1(i,1);
    end
end
% k1 = 0.5;

%%%%% checked with slides 56 and 59

% Initialize values for integration loop.
% theta = deg2rad(0:1:100);
r = zeros(10,length(gamma1));
thetadot = zeros(10,length(gamma1));
c = zeros(10,length(gamma1));
s = zeros(10,length(gamma1));
gamma = zeros(10,length(gamma1));
a_N = zeros(10,length(gamma1));
a_T = zeros(10,length(gamma1));
TOF_mat = zeros(length(gamma1),1);
theta = zeros(10,length(gamma1));

% The outer loop will analyze for every initial flight path angle.
for i = 1:length(gamma1)
    % Define some initial values for inner loop
    j = 1;
    theta(j,i) = theta0;
    % phi and k0 can then be determined
    
    % TEEEEEST!
    phi = acos(tan(gamma1(i,1))./(k1(i,1).*k2));
%     phi = -pi/2;
    
    % TEEEEEST!
    k0 = r1/(exp( k1(i,1).*sin(phi) ));
%     k0 = 150e6;

    c(j,i) = cos(k2*theta(j,i) + phi);
    s(j,i) = sin(k2*theta(j,i) + phi);
    r(j,i) = k0*exp(k1(i,1)*s(j,i));
    gamma(j,i) = gamma1(i,1);

    
    % V1.2 Acceleration at gamma(j,i)
    a_N(j,i) = ((n)*tan(gamma(j,i))/(2*cos(gamma(j,i)))) * ( (1/(tan(gamma(j,i))^2 + k1(i,1)*k2^2*s(j,i) + 1)) - ( (k2^2*(1-2*k1(i,1)*s(j,i)))/ (tan(gamma(j,i))^2 + k1(i,1)*k2^2*s(j,i) + 1)^2 ) );
    a_T(j,i) = a_N(j,i)*(mu/r(j,i)^2);
    
    % Now the determination of thetadot, a new theta, and a new gamma.
    thetadot(j,i) = sqrt( (mu/r(j,i).^3) / (tan(gamma(j,i))^2 + k1(i,1)*k2^2*s(j,i) + 1) );
    theta(j+1,i) = theta(j+1,i) + thetadot(j,i)*timestep*daylength;

    
    % The inner loop will integrate to obtain time of flight    
    while theta(j,i) < thetabar && j*timestep < TOF_req %  <== MODIFIED to be limited by thetabar instead of r2!
        j = j+1;

        c(j,i) = cos(k2*theta(j,i) + phi);
        s(j,i) = sin(k2*theta(j,i) + phi);
        r(j,i) = k0*exp(k1(i,1)*s(j,i));
        gamma(j,i) = atan(k1(i,1)*k2*c(j,i));
        
        % V1.2 Acceleration at gamma(j,i)
        a_N(j,i) = (n*tan(gamma(j,i))/(2*cos(gamma(j,i)))) * ( (1/(tan(gamma(j,i))^2 + k1(i,1)*k2^2*s(j,i) + 1)) - ( (k2^2*(1-2*k1(i,1)*s(j,i)))/ (tan(gamma(j,i))^2 + k1(i,1)*k2^2*s(j,i) + 1)^2 ) );
        a_T(j,i) = a_N(j,i).*(mu/r(j,i)^2);

        % Now the determination of thetadot, a new theta, and a new gamma.
        thetadot(j,i) = sqrt( (mu/r(j,i).^3) / (tan(gamma(j,i))^2 + k1(i,1)*k2^2*s(j,i) + 1) );
        theta(j+1,i) = theta(j,i) + thetadot(j,i)*timestep*daylength;
    end
        
    %%%%% Accelerations checked with slide 42
    
    % Now the interim matrix thingy can be determined, which will later be
    % summed to obtain the TOF.
    dtheta = thetadot.*timestep;
    interim = sqrt( r(:,i).^3.*( tan(gamma(:,i)).^2 + k1(i,1)*k2^2*s(:,i) + 1) ./ mu ) .*dtheta(:,i);
    TOF = sum(interim);
    TOF_mat(i,1) = TOF;
end

dV_mat = sum(a_T.*timestep.*daylength);

% % plot stuff
% figure
% [row,col,radius] = find(r(:,1));
% [row,col,angle] = find(theta(:,1));
% polar(angle,radius,'-')


