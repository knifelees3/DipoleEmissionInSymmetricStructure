Main_Basic_RhoPhi;
% angle=[pi/3,pi/2,1]';
% kx_grid_1d=reshape(kx_grid,num_kx*num_ky,1);
% ky_grid_1d=reshape(ky_grid,num_kx*num_ky,1);
% kxy_1d=[kx_grid_1d,ky_grid_1d];
% %%
% Pattern1d=nPattern_Cal_single(angle,kxy_1d);

%%
angle=[160/180*pi,26/180*pi,248/160*pi]';
%kphi_2=[kphi,kphi];
krho_2=linspace(0+1e-5,1+1e-5,num_l*2)*axis;
Pattern1d=nPattern_Cal_single_krho(angle,krho_2');

[kRho_grid,Phi_grid]=meshgrid(krho_2,kphi); % This definition is not so good; 
kx_grid_2=kRho_grid.*cos(Phi_grid);
ky_grid_2=kRho_grid.*sin(Phi_grid);
%Pattern=reshape(Pattern1d,num_l,num_ky);
%%
 %figure(1);pcolor(ky_grid_2,kx_grid_2,Pattern1d);shading interp; colormap jet; colorbar;
% %%
figure(2);plot(krho_2,Pattern1d,'*')

