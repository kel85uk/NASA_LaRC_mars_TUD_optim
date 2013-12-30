function[k2_final,N_final,gamma1_final,TOF_final,dV_trans,m_expo]=LT_EMtransfer(r1,r2,mu,TOF_req,timestep,payload,T,Isp)
% function[k2_final,N_final,gamma1_final,TOF_final,dV_trans,m_expo]=LT_EMtransfer(r1,r2,mu,TOF_req,timestep,payload,T,Isp)
%
% This function will optimize a transfer from Earth to Mars in heliocentric
% space, using the exposins2 function.
%
% [k2_final,N_final,gamma1_final,TOF_final,dV_trans,m_expo]=LT_EMtransfer(earth_sun_distance,mars_sun_distance,mu_sun,5*365,1,250000,0.05,4500)
%
%
% Written by Harry Linskens, 16-May-2013, modified from exposin2
%--------------------------------------------------------------------------
format longg;

% Values are initalized
Psi = 180;              % [deg], random choice
n = 1;                  % acceleration


% INPUT
%--------------------------------------------------------------------------
% The variables N, k2, and gamma1 are put in vectors.
N_vec = 1:1:5;                    % [revs]
k2_vec = 0.01:0.01:0.15;                % [-]
gamma1 = deg2rad(0:0.1:1);  % [rad]
%--------------------------------------------------------------------------

tic
% Initialize the output matrices
TOF_matr = zeros(length(k2_vec),length(N_vec),length(gamma1));
dV_matr = zeros(length(k2_vec),length(N_vec),length(gamma1));
a_matr = zeros(Plength(k2_vec),length(N_vec),length(gamma1));
counter = 0;

% Initialize waitbar
h = waitbar(0,'Doing all kinds of fancy things...');
set(findobj(h,'type','patch'),'edgecolor',[0 0.7 0.8],'facecolor',[0 0.7 0.8])

% The optimizer will loop through all options
for i = 1:length(N_vec)
    for j = 1:length(k2_vec)
        for k = 1:length(gamma1)
            N = N_vec(1,i);             % [revs]
            k2 = k2_vec(1,j);           % [-]
            [~,~,TOF_mat,~,dV_mat,a_T] = exposin2(r1,r2,Psi,TOF_req,timestep,N,k2,mu,n,gamma1(1,k));
            
            % Store the relevant values
            TOF_matr(j,i,k) = TOF_mat;
            dV_matr(j,i,k) = dV_mat;
            a_matr(j,i,k) = max(a_T);
            
            % give counter
            counter = counter+1;
            progress = counter/(length(N_vec)*length(k2_vec)*length(gamma1));
            waitbar(progress,h)
        end
    end
end
runtime2 = toc
close(h)


%% V1.1 MASS ITERATION LOOP
% The actual thruster that will be used is the RIT-XT, developed from the
% RIT-22 for commercial space applications.
% Isp = 4500;         % [s]

% It is desirable to keep the thrust level under 50 mN, as this reduces the
% necessary beam power from roughly 4700 W to 2000 W! This is due to a
% maximum eclipse in the parking orbit of roughly 37.5 minutes and 69.4
% minutes in GEO (SMAD).
% T_2000 = 50e-3;         % [N]
% T_2400 = 75e-3;         % [N]
% T_4700 = 150e-3;        % [N]
% T = T_2000;

% Preliminary system mass. INPUT A dV AND ITERATE. Initial estimate of
% delta-V is given by difference in circular velocities
dV_est = sqrt(mu/r1) - sqrt(mu/r2);
m_dry = payload;                           % [kg]
m_expo = m_dry*exp(dV_est*1000/Isp/9.81);       % [kg]


% With these parameters selected, the non-feasible options can be
% eliminated. Also those with a time of flight of the limit or more,
% defined by a certain tolerance to account for machine errors.
tolerance = 0.001;

for i = 1:length(N_vec)
    for j = 1:length(k2_vec)
        for k = 1:length(gamma1)
            % Filter by time of flight
            if TOF_matr(j,i,k) >= TOF_req-tolerance
                TOF_matr(j,i,k) = NaN;
            end
            
            % check acceleration level
            if isnan(TOF_matr(j,i,k)) == 1
                a_matr(j,i,k) = NaN;
            elseif a_matr(j,i,k) > T/m_expo
                a_matr(j,i,k) = NaN;
            end
            
            % Eliminate dVs which can't be reached
            if isnan(a_matr(j,i,k)) == 1
                dV_matr(j,i,k) = NaN;
            end
        end
    end
end


%% FINDING FINAL LOW-THRUST TRAJECTORY
% Find the minimum delta-V, and the corresponding transfer parameters.
% V1.1: Post dV_trans processing MOVED TO AFTER LOOP
dV_trans = min(min(min(dV_matr)));


[I J] = find(dV_matr == dV_trans);

if J <= length(N_vec)
    index_N = J;
else
    index_N = J - ( floor(J/length(N_vec)) * length(N_vec) );
end
index_g = floor(J/length(N_vec)) + 1;

k2_final = k2_vec(I);                        
N_final = N_vec( index_N );                  % [revs]
gamma1_final = rad2deg(gamma1( index_g ));   % [rad]
TOF_final = TOF_matr(I,index_N,index_g);     % [days]


% With these transfer properties known, the acceleration profile can be
% constructed quickly, along with the trajectory.
[r,~,~,theta,~,a_T] = exposin2(r1,r2,Psi,TOF_final,timestep,N_final,k2_final,mu,n,gamma1_final);
t_vec = linspace(0,TOF_final,length(a_T));

% Plot the trajectory in the equatorial plane
figure
[row,col,radius] = find(r(:,1));
[row,col,angle] = find(theta(:,1));
polar(angle,radius,'b-')
grid on
hold on
theta_orb = deg2rad(0:1:360);
r_E = ones(1,length(theta_orb))*r1; 
r_M = ones(1,length(theta_orb))*r2;
polar(theta_orb,r_M,'r')
polar(theta_orb,r_E,'g')
hold off
legend('Final transfer trajectory','Mars orbit','Earth orbit')
xlabel('Radial distance [km]')

% Plot the acceleration as a function of time
figure
plot(t_vec,a_T)
hold on
plot(t_vec,ones(length(t_vec)).*(T/m_expo))
hold off
grid on

