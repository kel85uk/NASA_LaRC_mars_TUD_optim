% OTV tank sizing extra
%
% Considers the case when the tanks surround the BA330. Determines the
% necessary radius for each perimiter tank, and the corresponding volume.
%
% Written by Harry Linskens 30-May-2013
%--------------------------------------------------------------------------

N = 6:1:12;
alpha = deg2rad(360./N);
r1 = 3.5;                       % [m] radius of BA330 thing (7m diamater)
L = 12;                          % [m]

% Now, find the corresponding r2 vector
r2 = sin(alpha./2)*r1 ./ (1-sin(alpha./2));

% Compute the total volume of the tanks (spherical end caps)
V = ((4/3)*pi*r2.^3) + (pi*r2.^2.*(L-(2*r2)));      % [m3] single tank
V_tot = V.*N;                                       % [m3] total volume

r2
V
V_tot