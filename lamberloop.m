% luca's lambert loop

%Code uses basic mars paramters and premade Lambert targeter.

clear all
clc
close all
%Set grav parameter for 2 planets
mu1=398600;
mu2=42828;

%Set their radii in AU
r1=1;
r2=1.5234;

%Set starting and ending eccentricities
e1=0;
e2=0;

%Set starting and ending semi major axes
a1_park=450+6378.138;
a2_park=300+6378.137;

% %Set resolution
% lim1=20;
% lim2=20;

%Code goes through this range of angles and time intervals set here in
%degrees and days
Dtheta=50:5:180;
Dt=120:5:270;


h = waitbar(0,'Doing all kinds of fancy things...');
set(findobj(h,'type','patch'),'edgecolor',[0 0.7 0.8],'facecolor',[0 0.7 0.8])
counter=0;

for k1=1:length(Dtheta)
    for k2=1:length(Dt)
       
        dtheta=Dtheta(k1);
        dt=Dt(k2);
        
        [a,dVtot,dV0,dV3] = Lambert_targ(dtheta,dt,mu1,mu2,r1,r2,a1_park,a2_park,e1,e2);
        %[deltav1 deltav2 dVtot kep_state a e EE w r2]=lamb('m',dt,dtheta);
        
        
        DV(k1,k2)=dVtot;
        
        counter=counter+1;
        progress=counter/(length(Dtheta)*length(Dt));
        waitbar(progress,h);
    end
end
close(h)


DTH=meshgrid(Dtheta);
DT=meshgrid(Dt);

surf(Dt,Dtheta,DV)
xlabel('Travel time')
ylabel('Phase angle')
%%
[total_mass mass_ratio]=propmass(DV,450,55);

total_mass(total_mass>1000)=NaN;

t_halving=12;       %[g/m^2];
T_tresh=120;
t=radprotection(T_thres,DT,t_halving);




figure;
surf(Dt,Dtheta,total_mass)
xlabel('Travel time')
ylabel('Phase angle')



