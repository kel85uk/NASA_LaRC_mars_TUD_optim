function[radiation_shielding_mass,mmod_shielding_mass] = shielding(crew_size,T_space)
% function[radiation_shielding_mass,mmod_shielding_mass] = shielding(crew_size,T_space)



% MMOD Shielding
shielding_MMOD = 37.5;          % kg/m^2 source:DSE Nerine

% volume needed per crew member for 180 days
habitable_volume_per_crew = 22; % m^3,  source:Factors Impacting Habitable Volume Requirements Figure1: Averaged habitable volume curve

% sizing based on a cilinder of d*d
total_required_habitable_volume = crew_size * habitable_volume_per_crew; % m^3

%keep fixed diameter
d = 7;
l = total_required_habitable_volume * 4 /d^2;

% area [m^2]
area_habitat = 2*pi*d/2 *l + 2*(pi*(d/2)^2);

T_thresh = 120;                         % [days]
m_htl = 120;                            % [kg/m2] for lead
t_shield = m_htl*log2(T_space/T_thresh); 

mmod_shielding_mass = shielding_MMOD*area_habitat;

radiation_shielding_mass = t_shield*area_habitat;



