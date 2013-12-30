function [V1, V2] = lambert(mu,R1, R2, t, string)
%global mu
r1 = norm(R1);
r2 = norm(R2);
c12 = cross(R1, R2);
theta = acos(dot(R1,R2)/r1/r2);
if nargin < 4 || (~strcmp(string,'retro') & (~strcmp(string,'pro')))
string = 'pro';
fprintf('\n ** Prograde trajectory assumed.\n')
end
if strcmp(string,'pro')
if c12(3) <= 0
theta = 2*pi - theta;
end
elseif strcmp(string,'retro')
if c12(3) >= 0
theta = 2*pi - theta;
end
end
A = sin(theta)*sqrt(r1*r2/(1 - cos(theta)));
z = -100;
while F(z,t) < 0
z = z + 0.1;
end
tol = 1.e-8;
nmax = 5000;
ratio = 1;
n = 0;
while (abs(ratio) > tol) & (n <= nmax)
n = n + 1;
ratio = F(z,t)/dFdz(z);
z = z - ratio;
end
if n >= nmax
fprintf('\n\n **Number of iterations exceeds %g \n\n ',nmax)
end
f = 1 - y(z)/r1;
g = A*sqrt(y(z)/mu);
gdot = 1 - y(z)/r2;
V1 = 1/g*(R2 - f*R1);
V2 = 1/g*(gdot*R2 - R1);
return

function dum = y(z)
dum = r1 + r2 + A*(z*S(z) - 1)/sqrt(C(z));
end

function dum = F(z,t)
dum = (y(z)/C(z))^1.5*S(z) + A*sqrt(y(z)) - sqrt(mu)*t;
end

function dum = dFdz(z)
if z == 0
dum = sqrt(2)/40*y(0)^1.5 + A/8*(sqrt(y(0)) + A*sqrt(1/2/y(0)));
else
dum = (y(z)/C(z))^1.5*(1/2/z*(C(z) - 3*S(z)/2/C(z)) ...
+ 3*S(z)^2/4/C(z)) + A/8*(3*S(z)/C(z)*sqrt(y(z)) ...
+ A*sqrt(C(z)/y(z)));
end
end

function dum = C(z)
dum = stumpC(z);
end

function dum = S(z)
dum = stumpS(z);
end

function s = stumpS(z)
if z > 0
s = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
elseif z < 0
s = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
else
s = 1/6;
end
end

function c = stumpC(z)
if z > 0
c = (1 - cos(sqrt(z)))/z;
elseif z < 0
c = (cosh(sqrt(-z)) - 1)/(-z);
else
c = 1/2;
end
end

end %lambert
