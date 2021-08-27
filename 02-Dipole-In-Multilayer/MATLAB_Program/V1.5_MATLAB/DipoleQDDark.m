% Give the dipole moment for a give angle
function [p1,p2,p3]=DipoleQDDark(para)

alpha=para(1);
phi=para(2);
ratio=para(3);

% For Dipole 1
d1x = -cos(alpha) * cos(phi);
d1y = -cos(alpha) * sin(phi);
d1z =  sin(alpha);

% For dipole 2
d2x =  sin(phi);
d2y = -cos(phi);
d2z =  0;

% For dark axis
d3x=sin(alpha)*cos(phi);
d3y=sin(alpha)*sin(phi);
d3z=cos(alpha);

nor=sqrt(1^2+1^2+abs(ratio)^2);
p1 =[d1x, d1y, d1z]'/nor;
p2 = [d2x, d2y, d2z]'/nor;
p3 =ratio*[d3x, d3y, d3z]'/nor;
end