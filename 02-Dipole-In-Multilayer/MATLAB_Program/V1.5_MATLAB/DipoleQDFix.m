% Give the dipole moment for a give angle
% The dipole can be determined by two variables and we let the another
% dipole's angle be fixed.

function [p1,p2]=DipoleQDFix(alpha, phi)

% For Dipole 1
d1x = -cos(alpha) * cos(phi);
d1y = -cos(alpha) * sin(phi);
d1z = sin(alpha);
% For dipole 2
d2x =sin(phi);

d2y = -cos(phi);

d2z = 0;

p1 =[d1x, d1y, d1z]';
p2 = [d2x, d2y, d2z]';

end