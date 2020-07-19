%% To plot the simulated file
load('./Data/DifferenceReNor_RhoPhi19_1009_2345.mat')
[xx,yy,zz]=meshgrid(phi_1_mat,alpha_mat,phi_2_mat);

nDiff_Temp=(Diff_Rho).^(-1);
nDiff=(nDiff_Temp/max(max(max(nDiff_Temp)))); % larger nDiff, more match

% Name to show the plot content
char_name='Rho';
color_1=0.99;
color_2=0.95;
color_3=0.9;
color_4=0.8;
% The function plots
SeriesPlot;