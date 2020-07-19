% To get the theta and phi for a given dipole
function [theta,phi]=Cal_ThetaPhi(p)
theta=acos(p(3));
phi_1=asin(p(2)/sin(theta));
phi_2=acos(p(1)/sin(theta));
if phi_1>=0 && phi_2<=pi/2
    phi=phi_1;
elseif phi_1>=0 && phi_2>pi/2
        phi=phi_2;
elseif phi_1<0 && phi_2<=pi/2
        phi=2*pi+phi_1;
elseif phi_1<0 && phi_2>pi/2
        phi=pi-phi_1;
end
end