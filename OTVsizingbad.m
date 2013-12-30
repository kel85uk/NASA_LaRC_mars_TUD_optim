function[OTV_finalmass,m_tank,m_payload,m_prop] = OTVsizingbad(dV_req,days_to_mars,days_return,days_stay,crewsize,Isp)
% function[OTV_finalmass] = OTVsizing(dV_req,mission_time);
%
% This function will size the OTV. First, the dry components will be
% selected, after which the fuel will be iterated for tank mass.
%
% Written by Harry Linskens, 23-May-2013
%--------------------------------------------------------------------------

% First, the dry components
m_capsule = 7000;           % [kg] mass of a Dragon? capsule
m_hab = 20000;              % [kg] mass of a BA330? habitat
m_engine = 2470;            % [kg] mass of a J-2X? engine
no_engines = 2;
m_lander=14000;

% Next, the consumables for the crew will be calculated
[m_there,~,m_stay,~,m_back] = ECLSS_lifesupport(crewsize,days_to_mars,1,days_stay,1,days_return);
m_con = m_there + m_stay + m_back;


% Also, the radiation shielding will be computed
[rad mmod] = shielding(crewsize,days_to_mars);
m_shield = rad;


% Now, compute preliminary payload mass
m_payload = m_capsule + m_hab + (no_engines*m_engine) + m_con + m_shield + m_lander;


% With that done, compute the wet mass using Tsiolkovsky
[total_mass ~] = propmass(dV_req,Isp,m_payload);
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

mfrac_tank = 0.001; %(SIVB_tankage+SII_tankage)/2;      
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

OTV_finalmass = total_mass;
m_prop = OTV_finalmass - m_tank - m_payload;


OTV_finalmass = total_mass;
m_prop = OTV_finalmass - m_tank - m_payload;
%%%%%

