%% To see whether the data is well fitted
% the distribution in kx,ky, in rho , in phi will be shown

%% Theoretical calculations
Main_Basic;
Main_Green;
char_name='FitNLM';

% alpha0=A(1);
% phi_10=A(2);
% phi_20=A(3);
alpha0=160/180*pi;
phi_10=26/180*pi;
phi_20=248/180*pi;

alpha0_angle=alpha0/pi*180;
phi_10_angle=phi_10/pi*180;
phi_20_angle=phi_20/pi*180;
% 
[p1_0,p2_0]=DipoleQD(alpha0,phi_10,phi_20);
%p1_0=[1,1,0]';
%p2_0=[0,0,0]';

[theta_1,gamma_1]=Cal_ThetaPhi(p1_0);
[theta_2,gamma_2]=Cal_ThetaPhi(p2_0);

theta_1_angle=theta_1/pi*180;
gamma_1_angle=gamma_1/pi*180;
theta_2_angle=theta_2/pi*180;
gamma_2_angle=gamma_2/pi*180;

Pattern_Fit=Cal_Pattern_QD(num_kx,num_ky,p1_0,p2_0,GreenSUp,GreenPUp,thetaUp);

% Normalization and transform
% Normalize with sum
% nPattern_Fit=Pattern_Fit/sum(Pattern_Fit,'all');
% Nor_BFP_Cut=Nor_BFP_Cut/sum(Nor_BFP_Cut,'all');

% Normalize with maximum value
nPattern_Fit=Pattern_Fit/max(max(Pattern_Fit));
Nor_BFP_Cut=Nor_BFP_Cut/max(max(Nor_BFP_Cut));

[P_rho,P_phi]=Transform_RhoPhi_Interp(kx_grid,ky_grid,kx_grid_in,ky_grid_in,nPattern_Fit,krho);
[Exp_rho,Exp_phi]=Transform_RhoPhi_Interp(kx_grid,ky_grid,kx_grid_in,ky_grid_in,Nor_BFP_Cut,krho);

%% Plot in the same figure
figure(5)
suptitlename=['\alpha=',num2str(alpha0_angle),'^{o}  \phi_{1}=',num2str(phi_10_angle),'^{o}  \phi_{2}=',num2str(phi_20_angle),'^{o}'];
suptitle(num2str(suptitlename));
subplot(2,2,1)
pcolor(kx_grid/k0,ky_grid/k0,nPattern_Fit);shading interp;colormap jet;colorbar;
xlabel('kx/k0')
ylabel('ky/k0')
title('Theoretical Fit')

subplot(2,2,2)
pcolor(kx_grid/k0,ky_grid/k0,Nor_BFP_Cut);shading interp;colormap jet;colorbar;
title('Experimental Data')
xlabel('kx/k0')
ylabel('ky/k0')

subplot(2,2,3)
plot(krho,P_rho,'r-','markersize',9)
hold on
plot(krho,Exp_rho,'r*')
legend('Fit','Experiment','location','best')
title('Along \rho')
xlabel('k_{\rho}/k_{0}')
ylabel('Intensity (arb. units)')

subplot(2,2,4)
plot(kphi,P_phi,'r-','markersize',9)
hold on
plot(kphi,Exp_phi,'r*')
hold on
legend('Fit','Experiment','location','best')
title('Along \phi')
xlabel('\phi')
ylabel('Intensity (arb. units)')
set(gcf,'unit','centimeters','position',[5 5 28 20]);
title(['Slice plot ',num2str(char_name)])
%savefig(['./Figures/CompareFourFigs_',num2str(char_name),'_.fig']);