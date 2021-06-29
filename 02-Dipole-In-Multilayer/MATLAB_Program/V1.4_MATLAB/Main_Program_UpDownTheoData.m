% To obtain the theoretical data

% The Up Green Data
Main_Basic_Up;
Main_Green;
SUpGreen=GreenSUp;
PUpGreen=GreenPUp;
Uptheta=thetaUp;

% The Up Green Data
Main_Basic_Down;
Main_Green;
SDownGreen=GreenSUp;
PDownGreen=GreenPUp;
Downtheta=thetaUp;

%%  define the dipole moment
px=[1 0 0];
py=[0 1 0];
pz=[0 0 1];

%% Calculate the far field for dx,dy,zd dipole 
tic
UpEFardx=Cal_Field_Dipole(px,SUpGreen,PUpGreen,Uptheta);
DownEFardx=Cal_Field_Dipole(px,SDownGreen,PDownGreen,Downtheta);

UpEFardy=Cal_Field_Dipole(py,SUpGreen,PUpGreen,Uptheta);
DownEFardy=Cal_Field_Dipole(py,SDownGreen,PDownGreen,Downtheta);

UpEFardz=Cal_Field_Dipole(pz,SUpGreen,PUpGreen,Uptheta);
DownEFardz=Cal_Field_Dipole(pz,SDownGreen,PDownGreen,Downtheta);
toc
%% Calculate the pattern
patdxUp=(abs(UpEFardx(:,:,1)).^2+abs(UpEFardx(:,:,2)).^2+abs(UpEFardx(:,:,3)).^2)./abs(cos(Uptheta).^0);
patdxDown=(abs(DownEFardx(:,:,1)).^2+abs(DownEFardx(:,:,2)).^2+abs(DownEFardx(:,:,3)).^2)./abs(cos(Downtheta).^0);
patdyUp=(abs(UpEFardy(:,:,1)).^2+abs(UpEFardy(:,:,2)).^2+abs(UpEFardy(:,:,3)).^2)./abs(cos(Uptheta).^0);
patdyDown=(abs(DownEFardy(:,:,1)).^2+abs(DownEFardy(:,:,2)).^2+abs(DownEFardy(:,:,3)).^2)./abs(cos(Downtheta).^0);
patdzUp=(abs(UpEFardz(:,:,1)).^2+abs(UpEFardz(:,:,2)).^2+abs(UpEFardz(:,:,3)).^2)./abs(cos(Uptheta).^0);
patdzDown=(abs(DownEFardz(:,:,1)).^2+abs(DownEFardz(:,:,2)).^2+abs(DownEFardz(:,:,3)).^2)./abs(cos(Downtheta).^0);

%% Interpolate in rho and phi directions
numin=num_kx;
thetamatin=linspace(0,pi/2,numin);
phimatin=linspace(0,2*pi,numin);
[thetagridin,phigridin]=meshgrid(thetamatin,phimatin);

% ux_grid_in=sin(thetagridin).*cos(phigridin);
% uy_grid_in=sin(thetagridin).*sin(phigridin);
% urho_grid_in=sqrt(ux_grid_in.^2+ux_grid_in.^2);
% [patdxUprho,patdxUpphi]=Transform_RhoPhi_Interp(uxtheo_grid,uytheo_grid,ux_grid_in,uy_grid_in,patdxUp,urho_grid_in);
% [patdxDownrho,patdxDownphi]=Transform_RhoPhi_Interp(uxtheo_grid,uytheo_grid,ux_grid_in,uy_grid_in,patdxDown,urho_grid_in);
% [patdyUprho,patdyUpphi]=Transform_RhoPhi_Interp(uxtheo_grid,uytheo_grid,ux_grid_in,uy_grid_in,patdyUp,urho_grid_in);
% [patdyDownrho,patdyDownphi]=Transform_RhoPhi_Interp(uxtheo_grid,uytheo_grid,ux_grid_in,uy_grid_in,patdyDown,urho_grid_in);
% [patdzUprho,patdzUpphi]=Transform_RhoPhi_Interp(uxtheo_grid,uytheo_grid,ux_grid_in,uy_grid_in,patdzUp,urho_grid_in);
% [patdzDownrho,patdzDownphi]=Transform_RhoPhi_Interp(uxtheo_grid,uytheo_grid,ux_grid_in,uy_grid_in,patdzDown,urho_grid_in);
%%
figure(1)
subplot(231)
pcolor(uxtheo_grid,uytheo_grid,patdxUp);shading interp;colormap('jet');xlabel('ux');ylabel('uy');title('Dx Up');
subplot(232)
pcolor(uxtheo_grid,uytheo_grid,patdyUp);shading interp;colormap('jet');xlabel('ux');ylabel('uy');title('Dy Up');
subplot(233)
pcolor(uxtheo_grid,uytheo_grid,patdzUp);shading interp;colormap('jet');xlabel('ux');ylabel('uy');title('Dz Up');
subplot(234)
pcolor(uxtheo_grid,uytheo_grid,patdxDown);shading interp;colormap('jet');xlabel('ux');ylabel('uy');title('Dx Down');
subplot(235)
pcolor(uxtheo_grid,uytheo_grid,patdyDown);shading interp;colormap('jet');xlabel('ux');ylabel('uy');title('Dy Down');
subplot(236)
pcolor(uxtheo_grid,uytheo_grid,patdzDown);shading interp;colormap('jet');xlabel('ux');ylabel('uy');title('Dz Down');