
% Run the structure set and Green Tensor program
BasicStructureSet;
ShowStructure(DinLayer);
ShowStructure(DinLayerRever)

%%  define the dipole moment
px=[1 0 0];
py=[0 1 0];
pz=[0 0 1];

%% Calculate the far field for dx,dy,zd dipole 
% objective=1 means the pattern in objective
% objective=0 means the pattern in far field sphere surface
objective=0;
patdxUp=Cal_Pattern_1DDipole(px,DinLayer,objective);
patdyUp=Cal_Pattern_1DDipole(py,DinLayer,objective);
patdzUp=Cal_Pattern_1DDipole(pz,DinLayer,objective);
patdxDn=Cal_Pattern_1DDipole(px,DinLayerRever,objective);
patdyDn=Cal_Pattern_1DDipole(py,DinLayerRever,objective);
patdzDn=Cal_Pattern_1DDipole(pz,DinLayerRever,objective);


ux_grid=DinLayer.ux_grid;
uy_grid=DinLayer.uy_grid;

figure()
subplot(231)
pcolor(ux_grid,uy_grid,patdxUp);shading interp;colormap('jet');xlabel('ux');ylabel('uy');title('Dx Up');
subplot(232)
pcolor(ux_grid,uy_grid,patdyUp);shading interp;colormap('jet');xlabel('ux');ylabel('uy');title('Dy Up');
subplot(233)
pcolor(ux_grid,uy_grid,patdzUp);shading interp;colormap('jet');xlabel('ux');ylabel('uy');title('Dz Up');
subplot(234)
pcolor(ux_grid,uy_grid,patdxDn);shading interp;colormap('jet');xlabel('ux');ylabel('uy');title('Dx Down');
subplot(235)
pcolor(ux_grid,uy_grid,patdyDn);shading interp;colormap('jet');xlabel('ux');ylabel('uy');title('Dy Down');
subplot(236)
pcolor(ux_grid,uy_grid,patdzDn);shading interp;colormap('jet');xlabel('ux');ylabel('uy');title('Dz Down');
sgtitle('MATLAB Program')
