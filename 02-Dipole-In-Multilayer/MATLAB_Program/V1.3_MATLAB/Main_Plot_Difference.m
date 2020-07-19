%% To plot the simulated file
load('./Data/DifferenceReNor_19_1010_2043.mat')
[xx,yy,zz]=meshgrid(phi_1_mat,alpha_mat,phi_2_mat);

nDiff_Temp=(Difference).^(-1);
nDiff=(nDiff_Temp/max(max(max(nDiff_Temp)))); % larger nDiff, more match

% Name to show the plot content
char_name='XY';
color_1=0.99;
color_2=0.95;
color_3=0.9;
color_4=0.8;
% The function plots
SeriesPlot;