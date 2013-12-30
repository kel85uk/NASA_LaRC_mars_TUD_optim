function [state1 ]= cart2kep(state2)


x1=state2(1);
y1=state2(2);
z1=state2(3);
x1dot=state2(4);
y1dot=state2(5);
z1dot=state2(6);
format shortg

Me=5.972E24;
G=6.67E-11;
mu=3.98600441E14;
mu=1.327178E20;
Re=6378136;

R1=[x1;y1;z1;];                 %Position vector is composed
V1=[x1dot;y1dot;z1dot];         %Velocity vector is composed     

r1=norm(R1);                    %Position magnitude is computed
v1=norm(V1);                    %Velocity magnitude is computed

H1=cross(R1,V1);                %Angular momentum is found
N1=cross([0;0;1],H1);           %N vector is defined

a1=((2/r1)-((v1^2)/mu))^-1;             %Semi major axis calculated
E1=((1/mu)*cross(V1,H1))-(1/r1)*R1;     %Eccentricity vector found
e1=norm(E1);                            %Eccentricity magnitude found

i1=acos(H1(3)/norm(H1));                %Inclination calculated

Nx=N1(1);                               %X component of N identified
Ny=N1(2);                               %Y component of N identified

Nxy=sqrt(Nx^2+Ny^2);                    %Magnitude of N vector found


RAAN1=atan2(Ny/Nxy,Nx/Nxy);             %RAAN is calculated
if RAAN1<0                              %Ensure RAAN is between 0 and 360
    RAAN1=RAAN1+2*pi;
end

w1=acos(dot(E1*(1/norm(E1)),N1*(1/norm(N1))));  %Find argument of periapsis
                                                %assuming it is positive
if dot(cross((N1*(1/norm(N1))),E1),H1)<0             %Should the condition not be met
    w1=-w1;                           %Invert the sign
end


theta1=acos(dot(R1*(1/norm(R1)),E1*(1/norm(E1))));  %Find true anomaly
                                                    %Assuming it positive
if dot(cross(E1,R1),H1)<0               %Should the condition not be met
    theta1=-theta1;                     %inver the sign
end


state1=[a1;e1;i1;RAAN1;w1;theta1];     %Compose the state vector
%state1(3:6)=state1(3:6)*180./pi;       %Put angles in degrees

