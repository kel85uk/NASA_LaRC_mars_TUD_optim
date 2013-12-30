function [OTV_finalmass,m_tank,m_dry,m_prop,m_shield,acc_max,acc_min,output_add] = OTVsizing(dV_req1,dV_req2,days_to_mars,days_stay,days_return,crewsize,Isp)
% function[OTV_finalmass] = OTVsizing(dV_req,mission_time);
%
% This function will size the OTV. First, the dry components will be
% selected, after which the fuel will be iterated for tank mass.
%
% Inputs: dV in km/s, Isp in s
%
% Written by Harry Linskens, 23-May-2013
%--------------------------------------------------------------------------
%format shortg

% First, the dry components
m_capsule = 7000;           % [kg] mass of a Dragon? capsule
m_hab = 20000;              % [kg] mass of a BA330? habitat
m_engine = 2470;            % [kg] mass of a J-2X? engine
no_engines = 2;


% Next, the consumables for the crew will be calculated
[m_there,~,m_stay,~,m_back] = ECLSS_lifesupport(crewsize,days_to_mars,1,days_stay,1,days_return);
m_con = m_there + m_stay + m_back;


% Also, the radiation shielding will be computed
[rad mmod] = shielding(crewsize,days_to_mars+days_return);
% IF ASSUMED THAT ON THE WAY THERE PROPELLANT SHIELDS EVERYTHING
%[rad mmod] = shielding(crewsize,days_return);

% Also, assume that a fraction of the BA330 contributes to the radiation
% shielding
if rad >= (m_hab*0.5)
    m_shield = rad-(m_hab*0.5);
elseif rad < (m_hab*0.5)
    m_shield = 0;
end

%nn=input('Include radiation protection (1 = yes, 0 = no)? ');
nn = 1;
if nn==0
    m_shield=0;
end

% Now, compute preliminary payload mass
m_payload = m_capsule + m_hab + (no_engines*m_engine) + m_con + m_shield;

% From first leg
[total_mass1,m_tank1,m_prop1,acc_max1,acc_min1,output_add1] = propulsion_subs(m_payload,dV_req1,Isp,no_engines);

% From second leg
[total_mass2,m_tank2,m_prop2,acc_max2,acc_min2,output_add2] = propulsion_subs(m_payload,dV_req2,Isp,no_engines);

if (total_mass2 > total_mass1)
  f_v = [total_mass1,m_tank1,m_prop1,acc_max1,acc_min1];
  total_mass = f_v(1);
  m_tank = f_v(2);
  m_prop = f_v(3);
  acc_max = f_v(4);
  acc_min = f_v(5);
%  [total_mass,m_tank,m_prop,acc_max,acc_min] = f_v;
  output_add = output_add1;
else
  f_v = [total_mass2,m_tank2,m_prop2,acc_max2,acc_min2];
  total_mass = f_v(1);
  m_tank = f_v(2);
  m_prop = f_v(3);
  acc_max = f_v(4);
  acc_min = f_v(5);  
%  [total_mass,m_tank,m_prop,acc_max,acc_min] = f_v;
  output_add = output_add2;
end

OTV_finalmass = total_mass;
m_dry = total_mass - m_prop;
