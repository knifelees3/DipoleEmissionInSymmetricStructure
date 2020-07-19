%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$User Defined Parameters$$$$$$$$$$$$$$$$$$$$$$$
% ****************************************************************************************
% Define the basic partameters
c = 3e8;
WL0 = 580e-9;
omega = 2 * pi / WL0 * c;

nUp = 1;
nDn = 1;

% num_layer means the number of layers and num_dl means the coordinate for
%different interface
num_layer = nUp + nDn + 1;

%The number of interfaces
num_dl = num_layer - 1;

%Initialize the coordinate for each layer
%The dipole should be better in the coordinate z=0
dl = zeros(num_dl, 1);

dis = 200e-9;


dl(2) = 600e-9;
dl(1) = 0e-9;

% The position of the dipole
POSD = dl(nDn+1) - dis;

% Initialize the permitivity
Eplist = zeros(num_layer,1);
% Eplist(5)=1.5^2;
% Eplsit(4)=-7.3734 + 1.1016 * 1i;
% Eplist(3) = 1.85^2;
% Eplist(2) = -7.3734 + 1.1016 * 1i;
% Eplist(1) = 1;

Eplist(3) = 1.78^2;
Eplist(2) = 1.5^2;
Eplist(1) = 1;

num_kx = 201;
num_ky = 201;

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
NA=1.78;
ThetaMax=asin(NA/sqrt(Eplist(num_layer)));
axis = k0 * NA;

%# The length to calculate the far field. However this is not needed if the
%range of k_rho is properly chosen
dUpFar = WL0 * 500;
dDnFar = -WL0 * 500;

% In the polar coordinate regular theta
num_Rho=num_kx;
num_Phi=num_kx;

Theta_mat=linspace(0+1e-5,ThetaMax+1e-5, num_Rho);
Phi_mat=linspace(0+1e-5,2*pi+1e-5,num_Phi);
[Theta_grid,Phi_grid]=meshgrid(Theta_mat,Phi_mat);
kx_grid=axis*sin(Theta_grid).*cos(Phi_grid);
ky_grid=axis*sin(Theta_grid).*sin(Phi_grid);
krho_grid = sqrt(kx_grid.^2 + ky_grid.^2);
krho_grid0=krho_grid/k0;


%In the polar coordinate, regular krho

% num_Rho=num_kx;
% num_Phi=num_kx;
% kRho_mat=linspace(0+1e-5,1+1e-5,num_Rho)*axis;
% Phi_mat=linspace(0+1e-5,2*pi+1e-5,num_Phi);
% [kRho_grid,Phi_grid]=meshgrid(kRho_mat,Phi_mat); % This definition is not so good; 
% kx_grid=kRho_grid.*cos(Phi_grid);
% ky_grid=kRho_grid.*sin(Phi_grid);
% krho_grid = sqrt(kx_grid.^2 + ky_grid.^2);
% krho_grid0=krho_grid/k0;


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
num_l=(1+num_kx)/2;
krho=linspace(0+1e-5,1+1e-5,num_l)*axis;
krho0=krho/k0;
kphi=linspace(0,2*pi,num_l);
% Distribution of the angle phi
kphi_grid=zeros(num_kx,num_kx);
for l=1:num_kx
    for m=1:num_ky
        if kx_grid(l,m)>=0 && ky_grid(l,m)>=0 
            kphi_grid(l,m)=asin(ky_grid(l,m)/krho_grid(l,m));
        elseif kx_grid(l,m)>=0 && ky_grid(l,m)<0
            kphi_grid(l,m)=2*pi+asin(ky_grid(l,m)/krho_grid(l,m)); 
        elseif kx_grid(l,m)<0 && ky_grid(l,m)<0
            kphi_grid(l,m)=pi-asin(ky_grid(l,m)/krho_grid(l,m)); 
        else
            kphi_grid(l,m)=pi-asin(ky_grid(l,m)/krho_grid(l,m)); 
        end
    end
end
%% Load the experimental data
% ***************************************************************************************
% load('./Data/Nor_BFP_Cut.mat')
% [Exp_rho,Exp_phi]=Transform_RhoPhi(num_kx,num_l,krho_grid0,kphi_grid,krho0,kphi,Nor_BFP_Cut);
% *******************************************
showtext=strcat(datestr(now,	'yyyy-mm-dd HH:MM:SS'),': Basic Parameters Have Been Prepared \n');
fprintf(showtext);