function nPattern=nPattern_Cal_single(angle,kxy)
global num_dl Eplist dl nUp nDn POSD dUpFar kl
alpha=angle(1);
phi_1=angle(2);
phi_2=angle(3);
kx=kxy(:,1);
ky=kxy(:,2);
[numk,~]=size(kx);
[p1,p2]=DipoleQD(alpha,phi_1,phi_2);
Pattern1d=zeros(numk,1);
for l=1:numk
kz=sqrt(kl.^2-kx(l)^2-ky(l)^2);
thetaUp=acos(kz(nUp+nDn+1)./kl(nUp+nDn+1));
 [GreenSUp,GreenPUp]=Cal_Green_Single(num_dl, Eplist, dl, kz, nUp, nDn, kx(l), ky(l),POSD,dUpFar);
Pattern=Cal_Pattern_QD_single(p1,p2,GreenSUp,GreenPUp,thetaUp);
Pattern1d(l)=Pattern;
end
% Pattern1d=reshape(Pattern,num_kx*num_ky,1);
nPattern=Pattern1d/sum(Pattern1d,'all');
end