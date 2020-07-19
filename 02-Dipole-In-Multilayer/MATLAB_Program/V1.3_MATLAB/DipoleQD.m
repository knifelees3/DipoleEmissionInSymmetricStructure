% Give the dipole moment for a give angle
function [p1,p2]=DipoleQD(alpha, phi_1, phi_2)

% For Dipole 1
d1x = sin(alpha) * cos(phi_1);
d1y = sin(alpha) * sin(phi_1);
d1z = cos(alpha);
% For dipole 2
d2x = -cos(phi_2) * cos(alpha) * cos(phi_1) + sin(phi_2) * sin(phi_1);

d2y = -cos(phi_2) * cos(alpha) * sin(phi_1) - sin(phi_2) * cos(phi_1);

d2z = cos(phi_2) * sin(alpha);

p1 =[d1x, d1y, d1z]';
p2 = [d2x, d2y, d2z]';

end