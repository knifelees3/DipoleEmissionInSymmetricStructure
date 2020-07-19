%% Experimental data analysis in the rho and phi direction for every single data
%Load the experimental data and plot it in the rho and phi direction
Main_Basic
Nor_BFP_Cut=Nor_BFP_Cut/sum(Nor_BFP_Cut,'all');
[Exp_rho,Exp_phi]=Transform_RhoPhi_Interp(kx_grid,ky_grid,kx_grid_in,ky_grid_in,Nor_BFP_Cut,krho);
nExp_rho=Exp_rho/max(Exp_rho);
nExp_phi=Exp_phi/max(Exp_phi);
% Exp_rho=sum(Exp_rho,1)/sum(Exp_rho,'all');
% Exp_phi=sum(Exp_phi,2)/sum(Exp_phi,'all');
%% Calculate the pattern for the horizontal dipole and vertical dipole
Main_Basic;
Main_Green;

%% Vertical Dipole
p1=[0,0,1]';
p2=[0,0,0];
PatternVED=Cal_Pattern_QD(num_kx,num_ky,p1,p2,GreenSUp,GreenPUp,thetaUp);
nPVED=PatternVED/sum(PatternVED,'all');
%nPVED=PatternVED/max(max(PatternVED));
%% Horizontal dipole
p1=[1,0,0]';
p2=[0,1,0]';
PatternHED=Cal_Pattern_QD(num_kx,num_ky,p1,p2,GreenSUp,GreenPUp,thetaUp);
nPHED=PatternHED/sum(PatternHED,'all');
%nPHED=PatternHED/max(max(PatternHED));
%% Distribution of data along rho and phi
[VED_rho,VED_phi]=Transform_RhoPhi_Interp(kx_grid,ky_grid,kx_grid_in,ky_grid_in,nPVED,krho);
[HED_rho,HED_phi]=Transform_RhoPhi_Interp(kx_grid,ky_grid,kx_grid_in,ky_grid_in,nPHED,krho);

%% Renormalzie the data in rho and phi
nVED_rho=VED_rho/max(VED_rho);
nHED_rho=HED_rho/max(HED_rho);

nVED_phi=VED_phi/max(VED_phi);
nHED_phi=HED_phi/max(HED_phi);
%% Plot the data
%% Plot the 2D data
% figure(1)
% set(gcf,'unit','centimeters','position',[1 1 25 25]);
% subplot(221)
% pcolor(kx_grid/k0,ky_grid/k0,Nor_BFP_Cut);shading interp;colormap jet; colorbar;xlabel('kx');ylabel('ky');title('Experimental Data')
% subplot(222)
% pcolor(kx_grid/k0,ky_grid/k0,nPVED);shading interp;colormap jet; colorbar;xlabel('kx');ylabel('ky');title('Vertical Dipole')
% subplot(223)
% pcolor(kx_grid/k0,ky_grid/k0,nPHED);shading interp;colormap jet; colorbar;xlabel('kx');ylabel('ky');title('Horizontal Dipole')


% %plot the data along rho
% y_temp=linspace(0,1e-4,10);
% x_temp=ones(10);
% subplot(224)
% 
figure(2)
plot(krho/k0,VED_rho,'-')
hold on
plot(krho/k0,HED_rho,'-')
hold on
plot(krho/k0,Exp_rho,'S')
legend('Vertical Dipole','Horizontal Dipole','Experiment','location','best')
xlabel('k_{\rho}/k_{0}')
ylabel('intensity')
title('Distribution Along k_{\rho}')
%set(gca,'linewidth',1.5,'fontsize',15)
%% Plot data along phi
% figure(3)
% plot(kphi,VED_phi,'-')
% hold on
% plot(kphi,HED_phi,'-')
% hold on
% plot(kphi,Exp_phi,'S')
% hold on

%% Compare with Weiwang Xu et al. Nano lett (2017)
theta_rho=asin(krho/ke)/pi*180;
figure(4)
plot(theta_rho,nVED_rho,'-')
hold on
plot(theta_rho,nHED_rho,'-')
hold on
plot(theta_rho,nExp_rho,'S')
legend('Vertical Dipole','Horizontal Dipole','Experiment','location','best')
xlabel('k_{\rho}/k_{0}')
ylabel('intensity')
grid
title('Distribution Along k_{\rho}')