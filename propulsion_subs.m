function [total_mass,m_tank,m_prop,acc_max,acc_min,output_add] = propulsion_subs(m_payload,dV_req,Isp,no_engines)
% Propulsion subsystem

% Compute the wet mass using Tsiolkovsky (Initialize iteration, M_tank = 0)
total_mass = propmass(dV_req,Isp,m_payload);
m_prop = total_mass - m_payload;

% Now, it is assumed that the tankage of the propellant consists of a
% certain percentage of the propellant mass.
% mfrac_tank = 0.08;          % based on Saturn V, roughly. 0.10 for Centaur

% Better estimate based on the S-II and S-IVB, this time taking the mass of
% the J2 engines on these stages into account!
% http://history.nasa.gov/SP-4029/Apollo_18-19_Ground_Ignition_Weights.htm
SIVB_prop = 87200 + 18000;
SIVB_dry = 9559;
J2_mass = 1788;

SIVB_tankage = (SIVB_dry - J2_mass)/SIVB_prop;
SII_tankage = ((481000*0.076) - (5*J2_mass)) / (481000*0.924);

mfrac_tank = (SIVB_tankage+SII_tankage)/2;      
m_propextra = 10;
m_tank = m_prop*mfrac_tank;
m_tank2 = m_tank;

% Iteration loop
while abs(m_propextra) > 1
    m_dry = m_payload+m_tank;
    total_mass2 = propmass(dV_req,Isp,m_dry);
    
    % Compute the difference
    m_propextra = total_mass2 - total_mass - m_tank2;
    total_mass = total_mass2;
    m_tank2 = m_propextra*mfrac_tank;                % ADDITIONAL TANK MASS!
    m_tank = m_tank + m_tank2;
end
if (isnan(total_mass)||isinf(total_mass))
  total_mass = 1e10;
  m_tank = 0.05*total_mass;
end

m_prop = total_mass - (m_tank + m_payload);


%% ADDITIONAL stuff

% Determine maximum acceleration
T_engine = 1307000;         % [N] thrust of J-2X? engine
acc_min = (no_engines*T_engine)/(total_mass)/9.81;           % [g]
acc_max = (no_engines*T_engine)/(m_tank + m_payload)/9.81;      % [g]


% Determine size of propellant tanks
OF_engine = 5;              % [-] between 4.5-5.5 for J-2X
rho_LOX = 1141;             % [kg/m3]
rho_LH2 = 70.85;            % [kg/m3]

m_LOX = m_prop/(OF_engine+1)*OF_engine;             % [kg]
m_LH2 = m_prop/(OF_engine+1);                       % [kg]

V_LOX = m_LOX/rho_LOX;      % [m3]
V_LH2 = m_LH2/rho_LH2;      % [m3]

% Required number of stages for trip 1
n_stages1 = 1;

% Required number of stages for trip 2
n_stages2 = 1;

output_add = [m_LOX,V_LOX;m_LH2,V_LH2;n_stages1,n_stages2];

end

function total_mass=propmass(delta_V,Isp,mass_payload)

%Simply Tsiokolvsky formula, inputs required: %Delta V, Isp of rocket, Mass
%of payload to be transported. There are two outputs: mass ratio, and total
%initial mass of rocket. If only the first two inputs are chosen, only the
%mass ratio is output. In only the first input is given, an Isp of 380 is
%used. Units of delta V is [km/s].

g0=9.81;
delta_V=delta_V*1000;

mass_ratio=exp(delta_V/(Isp*g0));
    
total_mass=mass_ratio*mass_payload;

end
