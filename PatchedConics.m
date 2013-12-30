function[dV1_maneuver,dV2_maneuver,dV_tot,TOF,Vhyp2] = PatchedConics(mu1,mu2,r1,r2,a1_park,a2_park,e1_park,e2_park)
% function[dV_tot] = PatchedConics(mu1,mu2,r1,r2,a1_park,a2_park,e1_park,e2_park)
%
% This function will determine a patched conics interplanetary transfer
% technique, utilizing the Hohmann trajectory. Assumptions are that the
% planetary orbits are circular and coplanar, and that they orbit the SUN.
%
% mu1 = gravitational parameter of departure planet [km3/s2]
% mu2 = gravitational parameter of arrival planet [km3/s2]
% r1 = orbital radius of departure planet [AU]
% r2 = orbital radius of arrival planet [AU]
% a1_park = semi-major axis of departure planet parking orbit [km]
% a2_park = semi-major axis of arrival planet parking orbit [km]
% e1_park = eccentricity of of departure planet parking orbit
% e2_park = eccentricity of of arrival planet parking orbit
%
% [dV1_maneuver,dV2_maneuver,dV_tot,TOF] = PatchedConics(398600.441,42832,1.00000011,1.52366231,6378.136+400,5500,0,0.3)
%
%
% Written by Harry Linskens
%
% VERSION HISTORY
% V1.0: Written 28-03-2013
%--------------------------------------------------------------------------
close all; clc;
format longg;

%% HELIOCENTRIC SCALE
AU = 149.597871e6;          % [km]
mu_Sun = 1.3271243e11;      % [km3/s2]
r1 = r1*AU;                 % [km]
r2 = r2*AU;                 % [km]


% Use the Hohmann function to determine the necessary velocities.
[Vinf1,Vinf2,~,TOF] = Hohmann(mu_Sun,r1,r2,0,0);


%% PLANETOCENTRIC SCALE - DEPARTURE
% Assume the necessary burns occur at perigee. Then, compute escape
% velocity and hyperbolic perigee velocity.
Vesc1 = sqrt(2*mu1/(a1_park*(1-e1_park)));                      % [km/s]
Vper1 = sqrt(mu1* (2/(a1_park*(1-e1_park)) - 1/a1_park) );      % [km/s]
Vhyp1 = sqrt(Vesc1^2 + Vinf1^2);                                % [km/s]

% Now, the first maneuver delta-V can be determined
dV1_maneuver = Vhyp1 - Vper1;                                   % [km/s]


%% PLANETOCENTRIC SCALE - ARRIVAL
% Assume the necessary burns occur at perigee. Then, compute escape
% velocity and hyperbolic perigee velocity.
Vesc2 = sqrt(2*mu2/(a2_park*(1-e2_park)));                      % [km/s]
Vper2 = sqrt(mu2* (2/(a2_park*(1-e2_park)) - 1/a2_park) );      % [km/s]
Vhyp2 = sqrt(Vesc2^2 + Vinf2^2);                                % [km/s]

% Now, the final maneuver delta-V can be determined
dV2_maneuver = Vhyp2 - Vper2;                                   % [km/s]


%% FINAL PROCESSING
dV_tot = dV1_maneuver + dV2_maneuver;                           % [km/s]



