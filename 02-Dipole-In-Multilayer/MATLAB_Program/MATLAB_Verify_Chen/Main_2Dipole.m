Main_Basic;
Main_Green;
%%
%[p1,p2]=DipoleQD(pi/2,0,0);
p1=[0,0,1]';
p2=[0,0,0]';
Pattern1=Cal_Pattern_QD(num_kx,num_ky,p1,p2,GreenSUp,GreenPUp,thetaUp);
nPattern1=Pattern1/sum(Pattern1,'all');
figure(1);pcolor(kx_grid/k0,ky_grid/k0,nPatternTest);shading interp;colormap jet; colorbar;
%figure(2);pcolor(kx_grid/k0,ky_grid/k0,Nor_BFP_Cut/sum(Nor_BFP_Cut,'all'));shading interp;colormap jet;
%%
nPatternTest=PatternTest/max(max(PatternTest));
P_rho_rhophi=sum(nPatternTest,1).*sin(Theta_mat);
%P_rho_rhophi=nPatternTest(1,:);
P_phi_rhophi=sum(nPatternTest,2);
figure(2)
plot((Theta_mat)/pi*180,P_rho_rhophi/max(P_rho_rhophi),'r')
xlabel('Angle')
xlim([0 90])
grid
% figure(2)
