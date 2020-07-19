% This is to fit the pattern function with the experimental data
Main_Basic;

kx_grid_1d=reshape(kx_grid,num_kx*num_ky,1);
ky_grid_1d=reshape(ky_grid,num_kx*num_ky,1);
kxy_1d=[kx_grid_1d,ky_grid_1d];

Nor_BFP_Cut_1d=reshape(Nor_BFP_Cut,num_kx*num_ky,1)/sum(Nor_BFP_Cut,'all');

modelfun=@nPattern_Cal_single;

angle0=[90/180*pi,90/180*pi,90/180*pi]';

opts = statset('Display','iter','TolFun',1e-20);
mdl=fitnlm(kxy_1d,Nor_BFP_Cut_1d,modelfun,angle0,'options',opts);

% fetch the fitted angle
CoeAngle=mdl.Coefficients.Estimate;
%% Calculate an arbitrary theoretical case as a test
% Main_Basic;
% Main_Green;
% 
% alpha0=0;
% phi_10=0;
% phi_20=0;
% 
% alpha0_angle=alpha0/pi*180;
% phi_10_angle=phi_10/pi*180;
% phi_20_angle=phi_20/pi*180;
% 
% [p1_0,p2_0]=DipoleQD(alpha0,phi_10,phi_20);
% 
% 
% P_theo=Cal_Pattern_QD(num_kx,num_ky,p1_0,p2_0,GreenSUp,GreenPUp,thetaUp);
% 
% n_theo=P_theo/sum(P_theo,'all');
% 
% n_theo_1d=reshape(n_theo,num_kx*num_ky,1);

%%
% opts = statset('Display','iter','TolFun',1e-20);
% mdl=fitnlm(kxy_1d,n_theo_1d,modelfun,angle0,'options',opts);