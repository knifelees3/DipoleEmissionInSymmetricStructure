%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$User Defined Parameters$$$$$$$$$$$$$$$$$$$$$$$
% ****************************************************************************************

% Gloabal variable definition 
global num_dl Eplist dl nUp nDn POSD dUpFar kl krho kphi
% Define the basic partameters
c = 3e8;
WL0 = 532e-9;
omega = 2 * pi / WL0 * c;

nUp = 2;
nDn = 2;

% num_layer means the number of layers and num_dl means the coordinate for
%different interface
num_layer = nUp + nDn + 1;

%The number of interfaces
num_dl = num_layer - 1;

%Initialize the coordinate for each layer
%The dipole should be better in the coordinate z=0
dl = zeros(num_dl, 1);

dis = 5e-9;
dl(4) = 200e-9;
dl(3) = 100e-9;
dl(2) = -100e-9;
dl(1) = -200e-9;

% The position of the dipole
POSD = dl(nDn+1) - dis;

% Initialize the permitivity
Eplist = zeros(num_layer);

Eplist(5) = 1;
Eplist(4) = 1;
Eplist(3) = 1;
Eplist(2) = 1.518^2;
Eplist(1) = 1.518^2;


num_kx = 155;
num_ky = 155;

% Definition of the 2D dipole orientation
num_alpha=91;
num_phi_1=181;
num_phi_2=181;

alpha_mat=linspace(0,pi,num_alpha);
phi_1_mat=linspace(0,pi*2,num_phi_1);
phi_2_mat=linspace(0,pi*2,num_phi_2);

%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$Derived Parameters$$$$$$$$$$$$$$$$$$$$$$$$$$
% ***************************************************************************************
% The range of the wave vector
k0 = 2 * pi / WL0;
kl = k0 * sqrt(Eplist);
ke=k0*sqrt(Eplist(num_layer));
NA=1.35;
axis = k0 * NA;

%# The length to calculate the far field. However this is not needed if the
%range of k_rho is properly chosen
dUpFar = WL0 * 500;
dDnFar = -WL0 * 500;

% In the polar coordinate regular theta
% num_Rho=num_kx;
% num_Phi=num_kx;
% ThetaMax=asin(NA/Eplist(num_layer));
% Theta_mat=linspace(0+1e-5,ThetaMax+1e-5, num_Rho);
% Phi_mat=linspace(0+1e-5,2*pi+1e-5,num_Phi);
% [Theta_grid,Phi_grid]=meshgrid(Theta_mat,Phi_mat);
% kx_grid=axis*sin(Theta_grid).*cos(Phi_grid);
% ky_grid=axis*sin(Theta_grid).*sin(Phi_grid);
% krho_grid = sqrt(kx_grid.^2 + ky_grid.^2);
% krho_grid0=krho_grid/k0;

% In the polar coordinate, regular krho
num_Rho=num_kx;
num_Phi=num_kx;
kRho_mat=linspace(0+1e-5,1+1e-5,num_Rho)*axis;
Phi_mat=linspace(0+1e-5,2*pi+1e-5,num_Phi);
[kRho_grid,Phi_grid]=meshgrid(kRho_mat,Phi_mat); % This definition is not so good; 
kx_grid=kRho_grid.*cos(Phi_grid);
ky_grid=kRho_grid.*sin(Phi_grid);
krho_grid = sqrt(kx_grid.^2 + ky_grid.^2);
krho_grid0=krho_grid/k0;


% In the cartisian coordinate
% kx = linspace(-1.0+1e-5, 1.0+1e-5, num_kx) * axis;
% ky = linspace(-1.0+1e-5, 1.0+1e-5, num_ky) * axis;
% [kx_grid, ky_grid] = meshgrid(ky, kx);
% krho_grid = sqrt(kx_grid.^2 + ky_grid.^2);
% krho_grid0=krho_grid/k0;

% The Z component of the wave vector
klz = zeros(num_kx, num_ky, num_layer);
theta = zeros(num_kx, num_ky, num_layer);

for l =1:num_layer
    klz(:, :, l) = sqrt(kl(l).^2 - krho_grid.^2);
    theta(:, :, l) = asin(krho_grid./kl(l));
end
thetaUp=theta(:,:,num_layer);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To remap the data in the polar coordinate
% num_l=(num_kx+1)/2;
num_l=num_kx;
kRho_mat_in=linspace(0+1e-5,1+1e-5,num_l)*axis;
Phi_mat_in=linspace(0+1e-5,2*pi+1e-5,num_l);
krho=kRho_mat_in';
kphi=Phi_mat_in';
[kRho_grid_in,Phi_grid_in]=meshgrid(krho,kphi); % This definition is not so good; 
kx_grid_in=kRho_grid_in.*cos(Phi_grid_in);
ky_grid_in=kRho_grid_in.*sin(Phi_grid_in);

%% Load the experimental data
% ***************************************************************************************
load('./Data/Nor_BFP_Cut.mat')
% load('./Data/Exp_rho.mat')
% load('./Data/Exp_phi.mat')
% Nor_BFP_Cut=Nor_BFP_Cut/sum(Nor_BFP_Cut,'all');
% [Exp_rho,Exp_phi]=Transform_RhoPhi_Interp(kx_grid,ky_grid,kx_grid_in,ky_grid_in,Nor_BFP_Cut);
% *******************************************
showtext=strcat(datestr(now,	'yyyy-mm-dd HH:MM:SS'),': Basic Parameters Have Been Prepared \n');
fprintf(showtext);