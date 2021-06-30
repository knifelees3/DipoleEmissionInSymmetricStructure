% Basic parameters of the structures.
% This program only calculate the far field pattern in the up direction.
% To obatin the far field in the sub


% Gloabal variable definition 
global num_dl Eplist dl nUp nDn POSD dUpFar kl krho kphi dDnFar
% Define the basic partameters
c = 3e8;
WL0 = 655e-9;
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

dis = 17e-9;


dl(2) = 100e-9;
dl(1) = 0e-9;

% The position of the dipole
POSD = dl(nDn+1) - dis;

% Initialize the permitivity
Eplist = zeros(num_layer);


% Dipole in two interface

Eplist(3) = 1.4564^2;
Eplist(2) = 1;
Eplist(1) = 1;

num_kx = 200;
num_ky = 200;


%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$Derived Parameters$$$$$$$$$$$$$$$$$$$$$$$$$$
% ***************************************************************************************
% The range of the wave vector
k0 = 2 * pi / WL0;
kl = k0 * sqrt(Eplist);
ke=k0*sqrt(Eplist(num_layer));
NA=1.4579;
k0NA= k0 * NA;

%# The length to calculate the far field. However this is not needed if the
%range of k_rho is properly chosen
% dUpFar = WL0 * 500;
% dDnFar = -WL0 * 500;
dUpFar = 1;
dDnFar = -1;


% In the cartisian coordinate
kx = linspace(-1.0+1e-5, 1.0+1e-5, num_kx) * k0NA;
ky = linspace(-1.0+1e-5, 1.0+1e-5, num_ky) * k0NA;
[kx_grid, ky_grid] = meshgrid(ky, kx);
uxtheo_grid=kx_grid/k0NA;
uytheo_grid=ky_grid/k0NA;

% The Z component of the wave vector
klz = zeros(num_kx, num_ky, num_layer);
theta = zeros(num_kx, num_ky, num_layer);
krho_grid=sqrt(kx_grid.^2+ky_grid.^2);
for l =1:num_layer
    klz(:, :, l) = sqrt(kl(l).^2 - krho_grid.^2);
    theta(:, :, l) = asin(krho_grid./kl(l));
end

thetaUp=theta(:,:,num_layer);

showtext=strcat(datestr(now,	'yyyy-mm-dd HH:MM:SS'),': Basic Parameters Have Been Prepared \n');
fprintf(showtext);