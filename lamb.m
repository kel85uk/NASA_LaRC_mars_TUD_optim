function [deltav1 deltav2 deltavtot kep_state a e EE w r2]=lamb(p2,dt,dtheta)
%For p2, input desired planet


%Define constants
AU=149597871000;

mu=1.327178E20;
mu_e=3.98600441E14;
Re=6378136;


%Based on planet, choose different values for constants
if p2=='j'
    mu_j=1.2669E17;
    Rj=71492000;
    r2=5.204267;
    r1_park=Re+0.1*Re;
    r2_park=Rj+0.1*Rj;
end

if p2=='m'
    mu_j=4.2828E13;
    Rj=3396200;
    r2=1.523679;
    %r2=1.5;
    r1_park=Re+0.1*Re;
    r2_park=Rj+0.1*Rj;    
end

if p2=='n'
    mu_j=6.896E15;
    Rj=22300000;
    r2=30;
    h1_park=200E3;
    h2_park=1000E3;
    r1_park=(Re+h1_park);
    r2_park=(Rj+h2_park);
end
 
%Define radii, angles and time in SI
r1=1;
r1=r1*AU;
r2=r2*AU;

dtheta=deg2rad(dtheta);
dt=dt*86164.1004;


%% Compute first intermediate parameters
c=sqrt(r1^2+r2^2-2*r1*r2*cos(dtheta));

s=(r1+r2+c)/2;


T=sqrt(8*mu/s^3)*dt;

q=(sqrt(r1*r2)/s)*cos(dtheta/2);


%% Here an initial guess for x is found

x=0;

E_lam=x^2-1;
y=sqrt(abs(E_lam)); 
z=sqrt(1-q^2+q^2*x^2); 
f=y*(z-q*x);
g=x*z-q*E_lam;  

if E_lam<0
    d=atan2(f./sqrt(f.^2+g.^2),g./sqrt(f.^2+g.^2));
end

if E_lam>0
        d=log(f+g);
end

T0=2.*(x-q.*z-(d./y))./E_lam;

if T0>T;
    x0=T0*(T0-T)/(4*T);
end

if T0<T;
    x0=-(T-T0)/(T-T0+4);
end

if T0==T
    x0=0;
end

h=0.001;
x=x0;

%% The Newton method is run to converge to an x
for k=1:100
    
    
    E_lam=x^2-1;
    y=sqrt(abs(E_lam)); 
    z=sqrt(1-q^2+q^2*x^2); 
    f=y*(z-q*x);
    g=x*z-q*E_lam;  

    if E_lam<0
        d=atan2(f./sqrt(f.^2+g.^2),g./sqrt(f.^2+g.^2));
    end

    if E_lam>0
        d=log(f+g);
    end
    
    F=T-2.*(x-q.*z-(d./y))./E_lam;
    
    E_lam_h=(x+h)^2-1;
    y_h=sqrt(abs(E_lam_h)); 
    z_h=sqrt(1-q^2+q^2*(x+h)^2); 
    f_h=y_h*(z_h-q*(x+h));
    g_h=(x+h)*z_h-q*E_lam_h;  

    if E_lam_h<0
        d_h=atan2(f_h./sqrt(f_h.^2+g_h.^2),g_h./sqrt(f_h.^2+g_h.^2));
    end

    if E_lam_h>0
        d_h=log(f_h+g_h);
    end
    
    F_h=T-2.*((x+h)-q.*z_h-(d_h./y_h))./E_lam_h;
    
    E_lam_hh=(x-h)^2-1;
    y_hh=sqrt(abs(E_lam_hh)); 
    z_hh=sqrt(1-q^2+q^2*(x-h)^2); 
    f_hh=y_hh*(z_hh-q*(x-h));
    g_hh=(x-h)*z_hh-q*E_lam_hh;  

    if E_lam_hh<0
        d_hh=atan2(f_hh./sqrt(f_hh.^2+g_hh.^2),g_hh./sqrt(f_hh.^2+g_hh.^2));
    end

    if E_lam_hh>0
        d_hh=log(f_hh+g_hh);
    end
    
    F_hh=T-2.*((x-h)-q.*z_hh-(d_hh./y_hh))./E_lam_hh;

    Fd=(F_h-F_hh)/(2*h);
    
    
    x=x-F/Fd;
    xx(k)=x;
    
    
end

%% Find all orbital parameters

a=(s/2)/(1-x^2);


gamma=sqrt(mu*s/2);
rho=(r1-r2)/c;
sigma=2*sqrt(r1*r2/(c^2))*sin(dtheta/2);

z=sqrt(1-q.^2+q.^2.*x.^2);

v_r1=gamma*((q*z-x)-rho*(q*z+x))/r1;
v_r2=-gamma*((q*z-x)+rho*(q*z+x))/r2;
v_t1=gamma*sigma*(z+q*x)/r1;
v_t2=gamma*sigma*(z+q*x)/r2;


R1=[r1; 0; 0];

R2=[r2*cos(dtheta); r2*sin(dtheta); 0];

er1=R1./norm(R1);
er2=R2./norm(R2);

VR1=v_r1*er1;
VR2=v_r2*er2;

if mod(dtheta,2*pi)==pi
    
    n=[0;0;1];
else n=cross(er1,er2);
    
end

VT1=v_t1*cross(n,er1);
VT1=v_t1*[0; 1 ;0];


VT2=v_t2*cross(n,er2);
VT2=v_t2*[-sin(pi-dtheta); -cos(pi-dtheta);0];

Vsat1=VR1+VT1;
Vsat2=VR2+VT2;

v_e=sqrt(mu/r1);
v_j=sqrt(mu/r2);

V_e=[0; v_e; 0];
V_j=[-v_j*sin(dtheta); v_j*cos(dtheta); 0];


V_rel1=Vsat1-V_e;
V_rel2=Vsat2-V_j;

fpa1=atan(v_r1/v_t1);
fpa2=atan(v_r2/v_t2);

e=sqrt(1-(r1*norm(Vsat1)^2/mu)*(2-r1*norm(Vsat1)^2/mu)*cos(fpa1)^2);

p=a*(1-e^2);

state=[R1;Vsat1]';

kep_state=cart2kep(state);


%% Compute velocities during the orbital transfer


v_park1=sqrt(mu_e/r1_park);
v_park2=sqrt(mu_j/r2_park);

v1_inf=norm(V_rel1);
v1_hyp=sqrt(2*mu_e/r1_park+v1_inf^2);

deltav1=v1_hyp-v_park1;

v2_inf=norm(V_rel2);
v2_hyp=sqrt(2*mu_j/r2_park+v2_inf^2);
deltav2=v2_hyp-v_park2;
deltav2=v_park2-v2_hyp;

deltavtot=abs(deltav1)+abs(deltav2);


%% NOT NECESSARY FOR EXERCISE
%maybe variables to be used later on, angular momentum, eccentricity vector
%and a strange argument of periapsis 
%(see http://en.wikipedia.org/wiki/Argument_of_periapsis)

H=cross(R1,Vsat1);

EE=(1/mu)*cross(Vsat1,H)-(1/norm(R1))*R1;

w=acos(EE(1)/norm(EE));

