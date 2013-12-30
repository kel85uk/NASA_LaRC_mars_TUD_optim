function [R1 V1] = ephemeris_stuff(target1,observer1,e_t)

cspice_furnsh('/home/kklloh/Documents/KAIST/Mission_design/naif0008.tls');

cspice_furnsh('/home/kklloh/Documents/KAIST/Mission_design/pck00010.tpc');

cspice_furnsh('/home/kklloh/Documents/KAIST/Mission_design/jup230l.bsp');

cspice_furnsh('/home/kklloh/Documents/KAIST/Mission_design/ura095.bsp');

cspice_furnsh('/home/kklloh/Documents/KAIST/Mission_design/sat341.bsp');

cspice_furnsh('/home/kklloh/Documents/KAIST/Mission_design/nep085.bsp');

cspice_furnsh('/home/kklloh/Documents/KAIST/Mission_design/plu021.bsp');

cspice_furnsh('/home/kklloh/Documents/KAIST/Mission_design/DE405/de405.bsp');

AU = 149597870.691; %Definition from JPL http://neo.jpl.nasa.gov/glossary/au.html
mu_sun = 132712440017.987;
mu_earth = 398600.433;
mu_mars = 42828;
%mu_venus = 324858.599;
mu_venus = mu_mars;

format long; format compact;

frame    = 'ECLIPJ2000';
abcorrect = 'None';

e_t = ['JD ',num2str(e_t,30)];

e_t = cspice_str2et(e_t);;

%et1 = cspice_str2et(time1);
%et2 = cspice_str2et(time2);

starg1 = mice_spkezr(target1, e_t, frame, abcorrect, observer1);

%fprintf('The position of  : %s \n', target1);
%fprintf('As observed from : %s \n', observer1);
%fprintf('In the reference frame   : %s \n', frame);
%utc_epoch = cspice_et2utc(et1, 'C', 3);
%fprintf('At Launch (UTC)     : %s (', utc_epoch);
%fprintf('%s) \n',cspice_et2utc(et1,'J',10));

R1 = starg1.state(1:3); V1 = starg1.state(4:6);

%COE1 = coe_from_sv(R1,V1,mu_sun);

%starg2 = mice_spkezr(target2, et2, frame, abcorrect, observer2);

%fprintf('The position of  : %s \n', target2);
%fprintf('As observed from : %s \n', observer2);
%fprintf('In the reference frame   : %s \n', frame);
%utc_epoch = cspice_et2utc(et2, 'C', 3);
%fprintf('At Arrival (UTC)     : %s (', utc_epoch);
%fprintf('%s) \n',cspice_et2utc(et2,'J',10));

%R2 = starg2.state(1:3); V2 = starg2.state(4:6);

%COE2 = coe_from_svd(R2,V2,mu_sun);
% 
% X1 = [R1 V1]';
% X2 = [R2 V2]';

cspice_unload('/home/kklloh/Documents/KAIST/Mission_design/naif0008.tls');

cspice_unload('/home/kklloh/Documents/KAIST/Mission_design/pck00010.tpc');

cspice_unload('/home/kklloh/Documents/KAIST/Mission_design/jup230l.bsp');

cspice_unload('/home/kklloh/Documents/KAIST/Mission_design/ura095.bsp');

cspice_unload('/home/kklloh/Documents/KAIST/Mission_design/sat341.bsp');

cspice_unload('/home/kklloh/Documents/KAIST/Mission_design/nep085.bsp');

cspice_unload('/home/kklloh/Documents/KAIST/Mission_design/plu021.bsp');

cspice_unload('/home/kklloh/Documents/KAIST/Mission_design/DE405/de405.bsp');

return

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function coe = coe_from_svd(R,V,mu)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
vector of orbital elements [h e RA incl w TA a]
%}
% ---------------------------------------------
eps = 1.e-10;
r = norm(R);
v = norm(V);
vr = dot(R,V)/r;
H = cross(R,V);
h = norm(H);
%...Equation 4.7:
incl = acos(H(3)/h);
%...Equation 4.8:
N = cross([0 0 1],H);
n = norm(N);
%...Equation 4.9:
if n ~= 0
  RA = acos(N(1)/n);
  if N(2) < 0
    RA = 2*pi - RA;
  end
  else
    RA = 0;
end
%...Equation 4.10:
E = 1/mu*((v^2 - mu/r)*R - r*vr*V);
e = norm(E);
%...Equation 4.12 (incorporating the case e = 0):
if n ~= 0
  if e > eps
    w = acos(dot(N,E)/n/e);
    if E(3) < 0
      w = 2*pi - w;
    end
  else
    w = 0;
  end
  else
    w = 0;
end
%...Equation 4.13a (incorporating the case e = 0):
if e > eps
  TA = acos(dot(E,R)/e/r);
  if vr < 0
    TA = 2*pi - TA;
  end
  else
    cp = cross(N,R);
  if cp(3) >= 0
    TA = acos(dot(N,R)/n/r);
  else
    TA = 2*pi - acos(dot(N,R)/n/r);
  end
end
%...Equation 4.62 (a < 0 for a hyperbola):
a = h^2/mu/(1 - e^2);
coe = [h e RA*180/pi incl*180/pi w*180/pi TA*180/pi a];
end %coe_from_sv
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end
