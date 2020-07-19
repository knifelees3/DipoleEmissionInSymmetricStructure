%% Experimental data analysis in the rho and phi direction
%Load the experimental data and plot it in the rho and phi direction
% 2019 10 16 (data before 2019 10 16)
Main_Basic
load('./Data/Nor_BFP_Cut_many.mat')
num_data=20;
nExp_rho=zeros(num_kx,num_data);
nExp_phi=zeros(num_ky,num_data);

for l=1:num_data
    [Exp_rho,Exp_phi]=Transform_RhoPhi_Interp(kx_grid,ky_grid,kx_grid_in,ky_grid_in,Nor_BFP_Cut(:,:,l),krho);
    nExp_rho(:,l)=Exp_rho/max(Exp_rho);
    nExp_phi(:,l)=Exp_phi/max(Exp_phi);
end
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


%%
figure(4)
plot(krho/k0,nVED_rho,'-')
hold on
plot(krho/k0,nHED_rho,'-')
hold on
for l=1:num_data
    plot(krho/k0,nExp_rho(:,l),'*')
    hold on
end
legend('Vertical Dipole','Horizontal Dipole','location','best')
xlabel('k_{\rho}/k_{0}')
ylabel('intensity')
grid
title('Distribution Along k_{\rho}')

%% Plot single
figure(5)
plot(krho/k0,nVED_rho,'-')
hold on
plot(krho/k0,nHED_rho,'-')
hold on
plot(krho/k0,nExp_rho(:,4),'*')
hold on
legend('Vertical Dipole','Horizontal Dipole','Data Num (3)','location','best')
xlabel('k_{\rho}/k_{0}')
ylabel('intensity')
grid
title('Distribution Along k_{\rho}')