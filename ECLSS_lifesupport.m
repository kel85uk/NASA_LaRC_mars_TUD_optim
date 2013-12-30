function[required_IMLEO_with_ECLSS1,required_IMLEO_with_ECLSS2,required_IMLEO_with_ECLSS3,required_IMLEO_with_ECLSS4,required_IMLEO_with_ECLSS5] = ECLSS_lifesupport(crewsize,days_to_mars,days_descent,days_surface,days_ascent,days_return)
% function[bla] = ECLSS_lifesupport(crewsize,days_to_mars,days_descent,days_surface,days_ascent,days_return)
%
% % Lifesupport using Human Spaceflight HSMAD L.K. Pranke W.J. Larson, editor. Human Space
% % flight: Mission Analysis and Design. The
% % McGraw-Hill Companies, Inc
% 
% clear all, clc;
% 
% % input
% passengers = 6
% days = 365
% 
% %lifesupport item in kg (HSMAD)
% dryfood = 0.6 * passengers * days
% water = 2.91 * passengers * days
% oxygen = 0.2925 * passengers * days
% nitrogen = 0.5 * passengers * days
% 
% %lifesupport systems based on 6 persons HSMAD
% ECLSS_air = 3056 % [kg] 4 BMS, TCCS, OGA, 3x redundant
% ECLSS_water = 300 % [kg] VCD (Vapor Compression Distillation)  , 2x redundant
%  
%  % total ECLSS (HSMAD)
%  
%  total_ECLSS = dryfood + water + oxygen + nitrogen + ECLSS_air + ECLSS_water
%  
%  %%

%% 
%  % ECLSS calculation mass for a crew of 6
%  % Lifesupport using Human Mission to Mars (D.Rapp 2008) with recycling
%  % using table 4.3 page 90

 
%  % Input days
%  days_to_mars = 180; % number of days, index = 1 , transit to Mars
%  days_descent = 15;  % number of days, index = 2 , descent
%  days_surface = 14;  % number of days, index = 3 , stay on Mars
%  days_ascent = 15;   % number of days, index = 4 , ascent
%  days_return = 180;  % number of days, index = 5 , return to Earth
%  
%  % Input crew size
%  crewsize = 1;
 
 % Input gear ratio
 gear_ratio1 = 1.3;    % transit to Mars
 gear_ratio2 = 1.3;  % descent
 gear_ratio3 = 1.3;    % surface stay on Mars
 gear_ratio4 = 1.3;   % ascent
 gear_ratio5 = 1.3;   % return to Earth
 

 
 % Transit to Mars index = 1
 
 %WATER in kg
 water_requirement1 = 29/180 * days_to_mars * 1000 * crewsize/6; % kg (M_T)
 water_ECLSS_plant1 = 1.4 * 1000;
 water_ELCSS_recovery1 = 99; % percent
 water_ECLSS_backup_cache1 = (1 - water_ELCSS_recovery1/100) * water_requirement1;  % kg
 total_water_ECLSS1 = water_ECLSS_plant1+water_ECLSS_backup_cache1; % kg
 
 Mt_Mls_water_ratio1 = water_requirement1/total_water_ECLSS1;
 
 % AIR in kg
 air_requirement1 = 4/180 * days_to_mars * 1000 * crewsize/6;
 air_ECLSS_plant1 = 0.5 * 1000; % kg
 air_ECLSS_recovery1 = 83; % percent
 air_ECLSS_backup_cache1 = (1-air_ECLSS_recovery1/100)*air_requirement1;  %  kg
 total_air_ECLSS1 = air_ECLSS_plant1 + air_ECLSS_backup_cache1;
 
 Mt_Mls_air_ratio1 = air_requirement1/total_air_ECLSS1;
 
 %FOOD/WASTE in kg
 food1 = 1.6/180 * days_to_mars * 1000 * crewsize/6; % kg
 waste_disposal_materials1 = 0.5/180 * days_to_mars * 1000; % kg
 
 %TOTAL
 total_mass_delivered_to_mars1 = total_water_ECLSS1 + total_air_ECLSS1 + food1 + waste_disposal_materials1;
 required_IMLEO_with_ECLSS1 = total_mass_delivered_to_mars1 * gear_ratio1;
 
 
 %%%%%%%%%%%%%%%%%% DESCENT,  index = 2

 
 %DESCENT WATER in kg
 water_requirement2 = 2/15 * days_descent * 1000 * crewsize/6; % kg (M_T)
%  water_ECLSS_plant2 = 1.4 * 1000
%  water_ELCSS_recovery2 = 99 % percent
%  water_ECLSS_backup_cache2 = (1 - water_ELCSS_recovery2/100) * water_requirement2  % kg
 total_water_ECLSS2 = water_requirement2; % kg
 
 Mt_Mls_water_ratio2 = water_requirement2/total_water_ECLSS2;
 
 % DESCENT AIR in kg
 air_requirement2 = 0.9/15 * days_descent * 1000 * crewsize/6;
%  air_ECLSS_plant2 = 0.5 * 1000 % kg
%  air_ECLSS_recovery2 = 83 % percent
%  air_ECLSS_backup_cache2 = (1-air_ECLSS_recovery2/100)*air_requirement2  %  kg
 total_air_ECLSS2 = air_requirement2;
 
 Mt_Mls_air_ratio2 = air_requirement2/total_air_ECLSS2;
 
 % DESCENT FOOD/WASTE in kg
 food2 = 0.15/15 * days_descent * 1000 * crewsize/6; % kg
 waste_disposal_materials2 = 0.05/15 * days_descent * 1000 * crewsize/6; % kg
 
 % DESCENT TOTAL in kg
 total_mass_delivered_to_mars2 = total_water_ECLSS2 + total_air_ECLSS2 + food2 + waste_disposal_materials2;
 required_IMLEO_with_ECLSS2 = total_mass_delivered_to_mars2 * gear_ratio2;
 
 %%%%%%%%%%%%%%%%%%%%
 
 %SURFACE STAY index = 3
 

%SURFACE STAY WATER in kg
 water_requirement3 = 29/180 * days_surface * 1000 * crewsize/6; % kg (M_T)
 water_ECLSS_plant3 = 1.4 * 1000;
 water_ELCSS_recovery3 = 94; % percent
 water_ECLSS_backup_cache3 = (1 - water_ELCSS_recovery3/100) * water_requirement3;  % kg
 total_water_ECLSS3 = water_ECLSS_plant3+water_ECLSS_backup_cache3; % kg
 
 Mt_Mls_water_ratio3 = water_requirement3/total_water_ECLSS3;
 
 %SURFACE STAY AIR in kg
 air_requirement3 = 4/180 * days_surface * 1000 * crewsize/6;
 air_ECLSS_plant3 = 0.5 * 1000; % kg
 air_ECLSS_recovery3 = 76; % percent
 air_ECLSS_backup_cache3 = (1-air_ECLSS_recovery3/100)*air_requirement3;  %  kg
 total_air_ECLSS3 = air_ECLSS_plant3 + air_ECLSS_backup_cache3;
 
 Mt_Mls_air_ratio3 = air_requirement3/total_air_ECLSS3;
 
 %SURFACE STAY FOOD/WASTE in kg
 food3 = 1.6/180 * days_surface * 1000 * crewsize/6; % kg
 waste_disposal_materials3 = 0.5/180 * days_surface * 1000 * crewsize/6; % kg
 
 %SURFACE STAY TOTAL
 total_mass_delivered_to_mars3 = total_water_ECLSS3 + total_air_ECLSS3 + food3 + waste_disposal_materials3;
 required_IMLEO_with_ECLSS3 = total_mass_delivered_to_mars3 * gear_ratio3;

 %%%%%%%%%%%%%%%%%%%%
 % ASCENT,  index = 4

 
 %ASCENT WATER in kg
 water_requirement4 = 2/15 * days_ascent * 1000 * crewsize/6; % kg (M_T)
%  water_ECLSS_plant2 = 1.4 * 1000
%  water_ELCSS_recovery2 = 99 % percent
%  water_ECLSS_backup_cache2 = (1 - water_ELCSS_recovery2/100) * water_requirement2  % kg
 total_water_ECLSS4 = water_requirement4; % kg
 
 Mt_Mls_water_ratio4 = water_requirement4/total_water_ECLSS4;
 
 % ASCENT AIR in kg
 air_requirement4 = 0.9/15 * days_ascent * 1000 * crewsize/6;
%  air_ECLSS_plant2 = 0.5 * 1000 % kg
%  air_ECLSS_recovery2 = 83 % percent
%  air_ECLSS_backup_cache2 = (1-air_ECLSS_recovery2/100)*air_requirement2  %  kg
 total_air_ECLSS4 = air_requirement4;
 
 Mt_Mls_air_ratio4 = air_requirement4/total_air_ECLSS4;
 
 % ASCENT FOOD/WASTE in kg
 food4 = 0.15/15 * days_ascent * 1000 * crewsize/6; % kg
 waste_disposal_materials4 = 0.05/15 * days_ascent * 1000 * crewsize/6; % kg
 
 % ASCENT TOTAL in kg
 total_mass_delivered_to_mars4 = total_water_ECLSS4 + total_air_ECLSS4 + food4 + waste_disposal_materials4;
 required_IMLEO_with_ECLSS4 = total_mass_delivered_to_mars4 * gear_ratio4;
 
 %%%%%%%%%%%%%%%%%%%%
 % RETURN TO EARTH index = 5
 
 %RETURN WATER in kg
 water_requirement5 = 29/180 * days_return * 1000 * crewsize/6; % kg (M_T)
 water_ECLSS_plant5 = 1.4 * 1000;
 water_ELCSS_recovery5 = 99; % percent
 water_ECLSS_backup_cache5 = (1 - water_ELCSS_recovery5/100) * water_requirement5;  % kg
 total_water_ECLSS5 = water_ECLSS_plant5+water_ECLSS_backup_cache5; % kg
 
 Mt_Mls_water_ratio5 = water_requirement5/total_water_ECLSS5;
 
 %RETURN AIR in kg
 air_requirement5 = 4/180 * days_return * 1000 * crewsize/6;
 air_ECLSS_plant5 = 0.5 * 1000; % kg
 air_ECLSS_recovery5 = 83; % percent
 air_ECLSS_backup_cache5 = (1-air_ECLSS_recovery5/100)*air_requirement5;  %  kg
 total_air_ECLSS5 = air_ECLSS_plant5 + air_ECLSS_backup_cache5;
 
 Mt_Mls_air_ratio5 = air_requirement5/total_air_ECLSS5;
 
 %RETURN FOOD/WASTE in kg
 food5 = 1.6/180 * days_return * 1000 * crewsize/6; % kg
 waste_disposal_materials5 = 0.5/180 * days_return * 1000 * crewsize/6; % kg
 
 %RETURN TOTAL
 total_mass_delivered_to_mars5 = total_water_ECLSS5 + total_air_ECLSS5 + food5 + waste_disposal_materials5;
 required_IMLEO_with_ECLSS5 = total_mass_delivered_to_mars5 * gear_ratio5;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% TOTAL MISSION ECLSS MASS in kg 
total_mission_ECLSS_mass = required_IMLEO_with_ECLSS1 + required_IMLEO_with_ECLSS2 + required_IMLEO_with_ECLSS3 + required_IMLEO_with_ECLSS4 + required_IMLEO_with_ECLSS5;
 
data = [required_IMLEO_with_ECLSS1 , required_IMLEO_with_ECLSS2 , required_IMLEO_with_ECLSS3 , required_IMLEO_with_ECLSS4 , required_IMLEO_with_ECLSS5] ;
% figure;
% bar(data)
% title('Consumable mass budget per module')
% set(gca,'XTickLabel',{'Transit to Mars','Descent','Surface stay','Ascent','Return to Earth'})
% ylabel('mass in kg')
% grid on