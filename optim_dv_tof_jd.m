% Optimization study main program (Hybrid: GA + fmincon)
% Requires MATLAB Parallel & Global Optimization Toolboxes + NASA's ODTBX
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%	Some functions are written by members of the team (credited in their individual function files)
% kklloh

clear all; clc;

matlabpool open;
%global mu mu_1 mu_2 stay_return dt_stay DV Isp OTV_finalmass 
%global m_tank m_payload m_prop m_shield acc_max acc_min output_add

mu = 132712440017.987;
mu_1 = 398600.433;
mu_2 = 42828;
deg = pi/180;
au = 149597870.691;
RP1 = 6378.1; %km
RP2 = 3396;

RP1_SOI = 925000;
RP2_SOI = 577000;

stay_return = 1; % 1 - return trip with stay at P2, 0 - one way trip to P2
dt_stay = 10; % Number of days to stay in planet 2 inclusive of orbital transfer from P2 surface to P2 parking
Isp = 340; % Assuming LOX-LCH4 combination

%...Departure date init
year0 = 2030;
month0 = 1;
day0 = 1;
hour0 = 0;
minute0 = 0;
second0 = 1;

%...Departure date min
year1 = 2025;
month1 = 1;
day1 = 1;
hour1 = 0;
minute1 = 0;
second1 = 1;

%...Departure date max
year2 = 2040; 
month2 = 1; 
day2 = 1; 
hour2 = 0; 
minute2 = 0; 
second2 = 1;

jd_0 = J0(year0,month0,day0) + (hour0 + minute0/60 + second0/3600)/24;
jd_min = J0(year1,month1,day1) + (hour1 + minute1/60 + second1/3600)/24;
jd_max = J0(year2,month2,day2) + (hour2 + minute2/60 + second2/3600)/24;

% Flight times
dt1_min = 120.0; % Transfer days (P1 -> P2)
dt1_max = 175.0;
dt2_min = 120.0; % Transfer days (P2 -> P1)
dt2_max = 200.0;

h_p1 = 300.0; % Initial P1 parking orbit altitude
h_p1_min = 300.0;
h_p1_max = 300.0;
h_p2 = 450.0; % Initial P2 parking orbit altitude
h_p2_min = 450.0;
h_p2_max = 450.0;

ecc_p1 = 0.0; % Initial P1 parking orbit eccentricity
ecc_p1_min = 0.0;
ecc_p1_max = 0.0;
ecc_p2 = 0.0; % Initial P2 parking orbit eccentricity
ecc_p2_min = 0.0;
ecc_p2_max = 0.0;

X0 = [jd_0,175.0,200.0,h_p1,h_p2,ecc_p1,ecc_p2]; %Baseline design and initial design state
Xmax = [jd_max;dt1_max;dt2_max;h_p1_max;h_p2_max;ecc_p1_max;ecc_p2_max];
Xmin = [jd_min;dt1_min;dt2_min;h_p1_min;h_p2_min;ecc_p1_min;ecc_p2_min];

options = optimset('Algorithm','interior-point','Display','iter','FinDifftype','central','UseParallel','always','PlotFcns',{@optimplotx,@optimplotfval,@optimplotconstrviolation});

% Use GA to produce a likely global minimum solution for initial condition
sizepop = 100;
Xmax1 = Xmax';
Xmin1 = Xmin';
P0 = Xmin1(ones(sizepop,1),:) + (Xmax1(ones(sizepop,1),:) - Xmin1(ones(sizepop,1),:)).*rand(sizepop,length(X0));
optionsga = gaoptimset('Display','iter','PlotFcns',{@gaplotbestf,@gaplotexpectation},'CreationFcn',@gacreationlinearfeasible,'CrossoverFcn',@crossoverheuristic,'PopInitRange',[Xmin';Xmax'],'PopulationSize',sizepop,'InitialPopulation',P0,'UseParallel','always','TolFun',1,'Generations',50,'StallGenLimit', 20);
%optionsga = gaoptimset('Display','iter','PlotFcns',{@gaplotbestf,@gaplotexpectation},'PopInitRange',[Xmin';Xmax'],'TolFun',1,'Generations',50,'StallGenLimit',20);

[X,fval,exitflag,optim_output] = optim_ez(X0,optionsga,options,Xmin,Xmax,mu,mu_1,mu_2,stay_return,dt_stay,Isp);

[DV,DV1,DV2,delta_V_D1,delta_V_A1,delta_V_D2,delta_V_A2] = transfer_ellipse_optim(X,mu,mu_1,mu_2,stay_return,dt_stay,Isp); % Calculate again optimum delta-V and transfer times...

cspice_furnsh('/home/kklloh/Documents/KAIST/Mission_design/naif0008.tls');
e_t = ['JD ',num2str(X(1),30)];
e_t = cspice_str2et(e_t);
launch_day = cspice_timout(e_t,'YYYY-MM-DD HR:MN:SC.### ::RND');
cspice_unload('/home/kklloh/Documents/KAIST/Mission_design/naif0008.tls');

fprintf('Optimal launch conditions found for a %d-day stay return mission. \n',dt_stay);
fprintf('Launch time = %s \n',launch_day);
fprintf('Transfer time (Planet 1 -> Planet 2) = %f days \n',X(2));
fprintf('Transfer time (Planet 1 <- Planet 2) = %f days \n',X(3));
fprintf('Parking orbit periapsis altitude at Planet 1 = %f km \n',X(4));
fprintf('Parking orbit eccentricity at Planet 1 = %f \n',X(6));
fprintf('Parking orbit periapsis altitude at Planet 2 = %f km \n',X(5));
fprintf('Parking orbit eccentricity at Planet 2 = %f \n',X(7));
fprintf('Delta V [P1 park -> P1 SOI] (km/s) = %f \n',delta_V_D1);
fprintf('Delta V [P2 SOI -> P2 park] (km/s) = %f \n',delta_V_A1);
fprintf('Delta V [Planet 1 -> Planet 2] (km/s) = %f \n',DV1);
fprintf('Delta V [P2 SOI <- P2 park] (km/s) = %f \n',delta_V_D2);
fprintf('Delta V [P1 park <- P1 SOI] (km/s) = %f \n',delta_V_A2);
fprintf('Delta V [Planet 1 <- Planet 2] (km/s) = %f \n',DV2);
fprintf('Total delta V (km/s) = %f \n',DV);
%fprintf('OTV Specifications \n\n');
%fprintf('Total Mass (kg) = %f \n',OTV_finalmass);
%fprintf('Propellant Mass (kg) = %f \n',m_prop);
%fprintf('Tank Mass (kg) = %f \n',m_tank);
%fprintf('Shield Mass (kg) = %f \n',m_shield);

%[U] = meshfree_func_plot(X); % To plot the performance of the final design point found

matlabpool close
