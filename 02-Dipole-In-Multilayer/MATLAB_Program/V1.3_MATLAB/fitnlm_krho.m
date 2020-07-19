% This is to fit the pattern function with the experimental data
Main_Basic_RhoPhi;
Nor_Exp_rho=reshape(Exp_rho/sum(Exp_rho(1:num_kx-1),'all'),num_kx,1);

krho_1d=reshape(krho,num_kx,1);
tbl=table(krho_1d,Nor_Exp_rho);

modelfun=@nPattern_Cal_single_krho;

angle0=[90/180*pi,0/180*pi,120/160*pi]';

opts = statset('Display','iter','TolFun',1e-20);
mdl=fitnlm(krho_1d,Nor_Exp_rho,modelfun,angle0,'options',opts);

% fetch the fitted angle
CoeAngle=mdl.Coefficients.Estimate;
