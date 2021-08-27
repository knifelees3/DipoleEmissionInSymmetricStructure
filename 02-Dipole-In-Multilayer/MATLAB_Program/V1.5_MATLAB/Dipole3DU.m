function [d1,d2,d3]=Dipole3DU(para)
alpha=para(1);
beta=para(2);

theta1=para(3);
phi1=para(4);
phi2=para(5);
dx=[1;0;0];
dy=[0;alpha;0];
dz=[0;0;beta];
uvec=[sin(theta1)*cos(phi1);sin(theta1)*sin(phi1);cos(theta1)];


rotation=[[cos(phi2)+uvec(1)^2*(1-cos(phi2)),uvec(1)*uvec(2)*(1-cos(phi2))-uvec(3)*sin(phi2),uvec(1)*uvec(3)*(1-cos(phi2))+uvec(2)*sin(phi2)];...
    [uvec(2)*uvec(1)*(1-cos(phi2))+uvec(3)*sin(phi2),cos(phi2)+uvec(2)^2*(1-cos(phi2)),uvec(2)*uvec(3)*(1-cos(phi2))-uvec(1)*sin(phi2)];...
    [uvec(3)*uvec(1)*(1-cos(phi2))-uvec(2)*sin(phi2),uvec(3)*uvec(2)*(1-cos(phi2))+uvec(1)*sin(phi2),cos(phi2)+uvec(3)^2*(1-cos(phi2))]];
norm=sqrt(1^2+alpha^2+beta^2);
d1=rotation*dx/norm;
d2=rotation*dy/norm;
d3=rotation*dz/norm;
end