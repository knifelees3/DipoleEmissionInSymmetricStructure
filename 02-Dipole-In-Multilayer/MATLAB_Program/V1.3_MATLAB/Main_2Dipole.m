%% Test program: Plot the emission pattern of an arbitray orientated 2D dipole
Main_Basic_RhoPhi;
Main_Green;
%% definition of the orientation
%[p1,p2]=DipoleQD(pi/2,0,0);

% alpha0=90/180*pi;
% phi_10=90/180*pi;
% phi_20=90/180*pi;
% A(1)=0;
% A(2)=0;
% A(3)=0;
% alpha0=A(1);
% phi_10=A(2);
% phi_20=A(3);
% 
% alpha0_angle=alpha0/pi*180;
% phi_10_angle=phi_10/pi*180;
% phi_20_angle=phi_20/pi*180;
% % 
% [p1,p2]=DipoleQD(alpha0,phi_10,phi_20);
p1=[1,0,0];
p2=[1,0,0];

%% Calculate the pattern and normalize the data
PatternTest=Cal_Pattern_QD(num_kx,num_ky,p1,p2,GreenSUp,GreenPUp,thetaUp);
nPatternTest=PatternTest/sum(PatternTest,'all');
figure(1);pcolor(kx_grid/k0,ky_grid/k0,nPatternTest);shading interp;colormap jet; colorbar;
%figure(2);pcolor(kx_grid/k0,ky_grid/k0,Nor_BFP_Cut/sum(Nor_BFP_Cut,'all'));shading interp;colormap jet;
%% In polar coordinate, the distribution along rho and phi can be obtained 
% nPatternTest=PatternTest/max(max(PatternTest));
% P_rho_rhophi=sum(nPatternTest,1);
% P_phi_rhophi=sum(nPatternTest,2);
% figure(2)
% plot(P_rho_rhophi)
% figure(3)
% plot(P_phi_rhophi)
% %% Far field 3D pattern
% xx=-PatternTest.*sin(kRho_grid/k0).*cos(Phi_grid);
% yy=-PatternTest.*sin(kRho_grid/k0).*sin(Phi_grid);
% zz=-PatternTest.*cos(kRho_grid/k0);

%figure(1);surf(xx,yy,zz);shading interp;colormap jet;