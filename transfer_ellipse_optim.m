function [DV,DV1,DV2,delta_V_D1,delta_V_A1,delta_V_D2,delta_V_A2] = transfer_ellipse_optim(X_vec,mu,mu_1,mu_2,stay_return,dt_stay,Isp)
%global mu mu_1 mu_2 stay_return dt_stay DV Isp OTV_finalmass 
%global m_tank m_payload m_prop m_shield acc_max acc_min output_add
deg = pi/180;
au = 149597870.691;
RP1 = 6378.1; %km
RP2 = 3396;

RP1_SOI = 925000;
RP2_SOI = 577000;

%...Departure
planet_id1 = 'Earth';

%...Arrival 
planet_id2 = 'Mars'; 

jd_init = X_vec(1);
dt1 = X_vec(2);
dt2 = X_vec(3);
h_p1 = X_vec(4);
h_p2 = X_vec(5);
ecc_p1 = X_vec(6);
ecc_p2 = X_vec(7);

DV = 0;
DV1 = 0;
DV2 = 0;
V_inf_D_s = 0;
delta_V_D1 = 0;
delta_V_A1 = 0;
delta_V_D2 = 0;
delta_V_A2 = 0;

dt_req1 = dt1;
dt_req2 = dt2;
jd_init2 = jd_init + dt_req1;

[R1, V1] = ephemeris_stuff(planet_id1, 'Sun', jd_init);
[R2, V2] = ephemeris_stuff(planet_id2, 'Sun', jd_init2);


if (dt_req1 >= 0)
  [VD1 VA2] = lambert(mu,R1,R2,dt_req1*86400,'pro');
else
  [VD1 VA2] = lambert(mu,R1,R2,dt_req1*86400,'retro');
end

%coe_tf = coe_from_sv(R1,VD1,mu); % Transfer trajectory orbital elements

V_inf_D = VD1 - V1; % Spacecraft hyperbolic excess velocity upon exiting Planet 1 SOI
V_inf_D_s = norm(V_inf_D); % Spacecraft excess speed

V_inf_A = VA2 - V2; % Spacecraft hyperbolic excess velocity upon crossing Planet 2 SOI
V_inf_A_s = norm(V_inf_A); % Spacecraft excess speed

% Departure trajectory (Assume Hohmann from parking orbit to SOI exit)
rp_D = RP1 + h_p1; % Periapsis altitude of P1
ap_D = rp_D/(1-ecc_p1);
vp_D = sqrt(V_inf_D_s^2 + 2*mu_1/rp_D); % Speed of spacecraft at periapsis of the departure hyperbola
vc_D = sqrt(mu_1/rp_D*(1+ecc_p1)); % Speed of spacecraft in P1 parking orbit
T_p1 = ap_D^(2/3)*2*pi()/(sqrt(mu_1)); % Period of P1 parking orbit
delta_V_D1 = abs(vp_D - vc_D); % Delta-v requirement for departure

% Arrival trajectory (Assume Hohmann to parking orbit from SOI entry)
rp_A = RP2 + h_p2; % Periapsis altitude of P2
ap_A = rp_A/(1-ecc_p2);
vp_A = sqrt(V_inf_A_s^2 + 2*mu_2/rp_A);
vc_A = sqrt(mu_2/rp_A*(1+ecc_p2));
T_p2 = ap_A^(2/3)*2*pi()/(sqrt(mu_2)); % Period of P2 parking orbit
delta_V_A1 = abs(vp_A - vc_A); % Delta-v requirement for arrival

DV1 = delta_V_D1 + delta_V_A1;

%fprintf('Total delta V from planet 1 to 2 (km/s) = %f \n',DeltaV1);
%fprintf('Total time of flight (days) = %f \n',dt_req1);


% Departure from planet 2
DeltaV2 = 0;
et2 = jd_init2;
if (stay_return == 1)
  % Stay to do some stuff in planet 2
  jd2 = et2 + dt_stay;
  dt_tf_b = dt_req2;
  jd2_fin = jd2 + dt_tf_b;

  [R1, V1] = ephemeris_stuff(planet_id2, 'Sun', jd2);
  [R2, V2] = ephemeris_stuff(planet_id1, 'Sun', jd2_fin);

  if (dt_req2 >= 0)
    [VD1 VA2] = lambert(mu,R1,R2,dt_req2*86400,'pro');
  else
    [VD1 VA2] = lambert(mu,R1,R2,dt_req2*86400,'retro');
  end

%  coe_tf = coe_from_sv(R1,VD1,mu); % Transfer trajectory orbital elements

  V_inf_D = VD1 - V1; % Spacecraft hyperbolic excess velocity upon exiting Planet 1 SOI
  V_inf_D_s = norm(V_inf_D); % Spacecraft excess speed

  V_inf_A = VA2 - V2; % Spacecraft hyperbolic excess velocity upon crossing Planet 2 SOI
  V_inf_A_s = norm(V_inf_A); % Spacecraft excess speed

  % Departure trajectory (Assume Hohmann from parking orbit to SOI exit)
  rp_D = RP2 + h_p2; % Periapsis altitude of P2
  ap_D = rp_D/(1-ecc_p2);
  vp_D = sqrt(V_inf_D_s^2 + 2*mu_2/rp_D); % Speed of spacecraft at periapsis of the departure hyperbola
  vc_D = sqrt(mu_2/rp_D*(1+ecc_p2)); % Speed of spacecraft in P1 parking orbit
  T_p2 = ap_D^(2/3)*2*pi()/(sqrt(mu_2)); % Period of P1 parking orbit
  delta_V_D2 = abs(vp_D - vc_D); % Delta-v requirement for departure

  % Arrival trajectory (Assume Hohmann to parking orbit from SOI entry)
  rp_A = RP1 + h_p1; % Periapsis altitude of P1
  ap_A = rp_A/(1-ecc_p1);
  vp_A = sqrt(V_inf_A_s^2 + 2*mu_1/rp_A);
  vc_A = sqrt(mu_1/rp_A*(1+ecc_p1));
  T_p1 = ap_A^(2/3)*2*pi()/(sqrt(mu_1)); % Period of P2 parking orbit
  delta_V_A2 = abs(vp_A - vc_A); % Delta-v requirement for arrival

  DV2 = delta_V_D2 + delta_V_A2;
  %fprintf('Total delta V from planet 2 to 1 (km/s) = %f \n',DeltaV2);
  %fprintf('Total time of flight (days) = %f \n',dt_req2);

end
DeltaV = DV1 + DV2;
DV = DeltaV;
%fprintf('Total delta V (km/s) = %f \n',DeltaV);
%fprintf('Total time of flight (days) = %f \n',dt_req1 + dt_req2);

%DV = [DeltaV, DeltaV1, DeltaV2, V_inf_D_s^2, delta_V_D1, delta_V_A1, delta_V_D2, delta_V_A2];


function y = zero_to_360(x)
if x >= 360
x = x - fix(x/360)*360;
elseif x < 0
x = x - (fix(x/360) - 1)*360;
end
y = x;
end %zero_to_360

end
