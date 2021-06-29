function [d1,d2,d3]=Dipole3DSymbolic(para)
alpha=para(1);
beta=para(2);
phix=para(3);
phiy=para(4);
phiz=para(5);
dx=[1;0;0];
dy=[0;alpha;0];
dz=[0;0;beta];

norm=sqrt(1^2+alpha^2+beta^2);

rotation=[[cos(phiz)*cos(phiy),cos(phiz)*sin(phiy)*sin(phix)-sin(phiz)*cos(phix),cos(phiz)*sin(phiy)*cos(phix)+sin(phiz)*sin(phix)];...
    [sin(phiz)*cos(phiy),sin(phiz)*sin(phiy)*sin(phix)+cos(phiz)*cos(phix),sin(phiz)*sin(phiy)*cos(phix)-cos(phiz)*sin(phix)];...
    [-sin(phiy),cos(phiy)*sin(phix),cos(phiy)*cos(phix)]];

d1=rotation*dx/norm;
d2=rotation*dy/norm;
d3=rotation*dz/norm;
end