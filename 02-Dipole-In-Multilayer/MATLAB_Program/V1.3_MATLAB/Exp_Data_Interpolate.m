load('./Data/Nor_BFP_Cut.mat');

c = 3e8;
WL0 = 638e-9;
omega = 2 * pi / WL0 * c;

k0 = 2 * pi / WL0;
NA=1.4;
axis = k0 * NA;


num_kx = 155;
num_ky = 155;

kx = linspace(-1.0+1e-5, 1.0+1e-5, num_kx) * axis;
ky = linspace(-1.0+1e-5, 1.0+1e-5, num_ky) * axis;
[kx_grid, ky_grid] = meshgrid(ky, kx);
krho_grid = sqrt(kx_grid.^2 + ky_grid.^2);
krho_grid0=krho_grid/k0;
%%
num_Rho=num_kx;
num_Phi=num_kx;
kRho_mat=linspace(1+1e-5,1.4+1e-5,num_Rho)*axis/NA;
krho=kRho_mat/axis;
Phi_mat=linspace(0+1e-5,2*pi+1e-5,num_Phi);
[kRho_grid,Phi_grid]=meshgrid(kRho_mat,Phi_mat); % This definition is not so good; 
kx_grid_in=kRho_grid.*cos(Phi_grid);
ky_grid_in=kRho_grid.*sin(Phi_grid);

BFP_interp = interp2(kx_grid,ky_grid,Nor_BFP_Cut,kx_grid_in,ky_grid_in);

BFP_rho=sum(BFP_interp,1).*krho;
BFP_phi=sum(BFP_interp.*krho,2);
%% Handle the abnormal data
BFP_phi(78)=(BFP_phi(77)+BFP_phi(79))/2;
BFP_rho(155)=20;

%figure(1);pcolor(kx_grid/k0,ky_grid/k0,Nor_BFP_Cut);shading interp;colormap jet; colorbar;xlabel('kx');ylabel('ky');
%figure(1);pcolor(kx_grid_in/k0,ky_grid_in/k0,BFP_interp);shading interp;colormap jet; colorbar;xlabel('kx');ylabel('ky');


figure(2);plot(kRho_mat/k0,BFP_rho/max(BFP_rho),'r*');xlabel('k_{\rho}/k_{0}');
figure(3);plot(Phi_mat,BFP_phi/max(BFP_phi),'rs');xlabel('\phi');

%%
Exp_rho=BFP_rho/sum(BFP_rho(1:num_Rho-1),'all');
Exp_phi=BFP_phi/sum(BFP_phi(1:num_Rho),'all');
% save('./Data/Exp_rho.mat','Exp_rho');
% save('./Data/Exp_phi.mat','Exp_phi');