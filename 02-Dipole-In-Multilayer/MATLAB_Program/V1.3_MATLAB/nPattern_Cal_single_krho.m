function nPattern=nPattern_Cal_single_krho(angle,krho)
global num_dl Eplist dl nUp nDn POSD dUpFar kl kphi
alpha=angle(1);
phi_1=angle(2);
phi_2=angle(3);
[p1,p2]=DipoleQD(alpha,phi_1,phi_2);
[num_rho,~]=size(krho);
[num_phi,~]=size(kphi);


PatternRho=zeros(num_phi,num_rho);
for l=1:num_phi
    for m=1:num_rho
        kx=krho(m)*cos(kphi(l));
        ky=krho(m)*sin(kphi(l));
    kz=sqrt(kl.^2-kx^2-ky^2);
    thetaUp=acos(kz(nUp+nDn+1)./kl(nUp+nDn+1));
     [GreenSUp,GreenPUp]=Cal_Green_Single(num_dl, Eplist, dl, kz, nUp, nDn, kx, ky,POSD,dUpFar);
    PatternRho(l,m)=Cal_Pattern_QD_single(p1,p2,GreenSUp,GreenPUp,thetaUp);
    end
end
Pattern1d=reshape(sum(PatternRho,1),num_rho,1);
nPattern=Pattern1d/sum(Pattern1d,'all');
end