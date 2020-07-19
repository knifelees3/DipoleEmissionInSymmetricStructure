% This is to fit the pattern function with the experimental data
Main_Basic_RhoPhi;
Nor_Exp_phi=reshape(Exp_phi/sum(Exp_phi(1:num_kx-1),'all'),num_kx,1);

kphi_1d=reshape(kphi,num_kx,1);
tbl=table(kphi_1d,Nor_Exp_phi);

modelfun=@nPattern_Cal_single_kphi;

angle0=[160/180*pi,26/180*pi,248/160*pi]';

mdl=fitnlm(kphi_1d,Nor_Exp_phi,modelfun,angle0);

% fetch the fitted angle
CoeAngle=mdl.Coefficients.Estimate;
