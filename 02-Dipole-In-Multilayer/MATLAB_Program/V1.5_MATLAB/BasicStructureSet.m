% Basic parameters of the structures.
% This program only calculate the far field pattern in the up direction.
% To obatin the far field in the sub

nSub=1.4564;


% Gloabal variable definition 
% global num_dl Eplist dl nUp nDn POSD dUpFar kl krho kphi dDnFar
% User defined parameters
% ---------------------------------------------------------------------------------------
DinLayer.c = 3e8;
DinLayer.WL0 = 655e-9;

DinLayer.NA=nSub;

DinLayer.nUp = 1;
DinLayer.nDn = 1;

% num_layer means the number of layers and num_dl means the coordinate for
%different interface
DinLayer.num_layer = DinLayer.nUp + DinLayer.nDn + 1;

%The number of interfaces
DinLayer.num_dl = DinLayer.num_layer - 1;

%Initialize the coordinate for each layer
%The dipole should be better in the coordinate z=0
DinLayer.dl = zeros(DinLayer.num_dl, 1);

DinLayer.dis = 17e-9;


DinLayer.dl(2) = 100e-9;
DinLayer.dl(1) = 0e-9;

% The position of the dipole
DinLayer.POSD = DinLayer.dl(DinLayer.nDn) + DinLayer.dis;

% Initialize the permitivity
DinLayer.Eplist = zeros(DinLayer.num_layer,1);


% Dipole in two interface

DinLayer.Eplist(3) = 1;
DinLayer.Eplist(2) = 1;
DinLayer.Eplist(1) = nSub^2;

DinLayer.num_kx = 100;
DinLayer.num_ky = 100;







% Derived Parameters
% ---------------------------------------------------------------------------------------
% The range of the wave vector
DinLayer.k0 = 2 * pi / DinLayer.WL0;
DinLayer.kl = DinLayer.k0 * sqrt(DinLayer.Eplist);
DinLayer.ke=DinLayer.k0*sqrt(DinLayer.Eplist(DinLayer.num_layer));
DinLayer.k0NA= DinLayer.k0 * DinLayer.NA;

%# The length to calculate the far field. However this is not needed if the
%range of k_rho is properly chosen
% dUpFar = WL0 * 500;
% dDnFar = -WL0 * 500;
DinLayer.dUpFar = 1;
DinLayer.dDnFar = -1;


% In the cartisian coordinate
DinLayer.kx = linspace(-1.0+1e-5, 1.0+1e-5, DinLayer.num_kx) * DinLayer.k0NA;
DinLayer.ky = linspace(-1.0+1e-5, 1.0+1e-5, DinLayer.num_ky) * DinLayer.k0NA;
[DinLayer.kx_grid, DinLayer.ky_grid] = meshgrid(DinLayer.kx, DinLayer.ky);
DinLayer.ux_grid=DinLayer.kx_grid/DinLayer.k0NA;
DinLayer.uy_grid=DinLayer.ky_grid/DinLayer.k0NA;

% The Z component of the wave vector
DinLayer.klz = zeros(DinLayer.num_kx, DinLayer.num_ky, DinLayer.num_layer);
DinLayer.theta = zeros(DinLayer.num_kx, DinLayer.num_ky, DinLayer.num_layer);
DinLayer.krho_grid=sqrt(DinLayer.kx_grid.^2+DinLayer.ky_grid.^2);
for l =1:DinLayer.num_layer
    DinLayer.klz(:, :, l) = sqrt(DinLayer.kl(l).^2 - DinLayer.krho_grid.^2);
    DinLayer.theta(:, :, l) = asin(DinLayer.krho_grid./DinLayer.kl(l));
end

DinLayer.thetaUp=DinLayer.theta(:,:,DinLayer.num_layer);
DinLayerRever=StructureReSet(DinLayer);


showtext=strcat(datestr(now,    'yyyy-mm-dd HH:MM:SS'),': Basic Parameters Have Been Prepared \n');
fprintf(showtext);

DinLayer=Cal_Green_List(DinLayer);
DinLayerRever=Cal_Green_List(DinLayerRever);
