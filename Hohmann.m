function[dV1_final dV2_final dVtot TOF] = Hohmann(mu,r1,r2,i1,i2)
% function[dV1 dV2 dVtot] = Hohmann(mu,r1,r2,i1,i2)
%
% mu = gravitational parameter of dominant celestial body [km3/s2]
% r1 = initial orbit radius [km]
% r2 = final orbit radius [km]
% i1 = initial inclination [deg]
% i2 = final inclination [deg]
%
% This function will determine the associated delta-vs with the Hohmann
% transfer specified in the input. Negative apogee and perigee burns
% indicate decelerations! Total dV is output positive.
%
% Written by Harry Linskens, 14-March-2013
%--------------------------------------------------------------------------


% In-plane Hohmann burns
a_hoh = (r1+r2)/2;          % [km]
V1c = sqrt(mu/r1);                % [km/s]
V1h = sqrt(mu*(2/r1 - 1/a_hoh));  % [km/s]
V2c = sqrt(mu/r2);                % [km/s]
V2h = sqrt(mu*(2/r2 - 1/a_hoh));  % [km/s]
dV1 = V1h - V1c;                    % [km/s]
dV2 = V2c - V2h;                    % [km/s]
TOF = pi*sqrt(a_hoh^3/mu)/86400;    % [solar days]

% Loop to optimize plane change. First, the case is considered of a
% decrease in inclination.
if i1 > i2
    i_step = -0.1;               % [deg]    
    i_vec = deg2rad(i1:i_step:i2);
    dV_hoh = zeros(length(i_vec),3);
    
    for i = 1:length(i_vec)
        i_GTO = deg2rad(i1) - i_vec(i);
        dV1_tot = sqrt(V1c^2 + V1h^2 - 2*V1h*V1c*cos(i_vec(i)));
        dV2_tot = sqrt(V2c^2 + V2h^2 - 2*V2h*V2c*cos(i_GTO));
        dV_hoh(i,1) = dV1_tot;
        dV_hoh(i,2) = dV2_tot;
        dV_hoh(i,3) = dV1_tot + dV2_tot;
    end
    
    [dVtot index_GTO] = min(dV_hoh(:,3));
    dV1_final = dV_hoh(index_GTO,1);
    dV2_final = dV_hoh(index_GTO,2);
        
    
% Now, the case is considered that the inclination must be increased.
elseif i2 > i1
    i_step = 0.1;                % [deg]
    
    i_vec = deg2rad(i1:i_step:i2);
    dV_hoh = zeros(length(i_vec),3);
    
    for i = 1:length(i_vec)
        i_GTO = deg2rad(i1) - i_vec(i);
        dV1_tot = sqrt(V1c^2 + V1h^2 - 2*V1h*V1c*cos(i_vec(i)));
        dV2_tot = sqrt(V2c^2 + V2h^2 - 2*V2h*V2c*cos(i_GTO));
        dV_hoh(i,1) = dV1_tot;
        dV_hoh(i,2) = dV2_tot;
        dV_hoh(i,3) = dV1_tot + dV2_tot;
    end
    
    [dVtot index_GTO] = min(dV_hoh(:,3));
    dV1_final = dV_hoh(index_GTO,1);
    dV2_final = dV_hoh(index_GTO,2);
    
        
% Finally, if the inclinations are equal...
elseif i1 == i2
    dV1_final = dV1;            % [km/s]
    dV2_final = dV2;            % [km/s]
    dVtot = abs(dV1_final+dV2_final);    % [km/s]
end
