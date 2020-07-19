% num_alpha=91;
% num_phi_1=1;
% num_phi_2=1;
% 
% alpha_mat=linspace(0,pi,num_alpha);
% phi_1_mat=linspace(0,pi*2,num_phi_1);
% phi_2_mat=linspace(0,pi*2,num_phi_2);

Diff_Rho=zeros(num_alpha,num_phi_1,num_phi_2);
Diff_Phi=zeros(num_alpha,num_phi_1,num_phi_2);
% *******************************************
showtext=strcat(datestr(now,	'yyyy-mm-dd HH:MM:SS'),': Begin Parallel Compute \n');
fprintf(showtext);

% Normalize with the average of the data
Nor_BFP_Cut=Nor_BFP_Cut/sum(Nor_BFP_Cut,'all');

kernum=parpool(60);
parfor l=1:num_alpha
    showtext=strcat(datestr(now,'yyyy-mm-dd HH:MM:SS'),': current loop is ',num2str(l),'\n');
    fprintf(showtext);
    for m=1:num_phi_1
        for n=1:num_phi_2
            [p1,p2]=DipoleQD(alpha_mat(l),phi_1_mat(m),phi_2_mat(n));
            PatternMat=Cal_Pattern_QD(num_kx,num_ky,p1,p2,GreenSUp,GreenPUp,thetaUp);
            [P_rho,P_phi]=Transform_RhoPhi(num_kx,num_l,krho_grid0,kphi_grid,krho0,kphi,PatternMat);
            [Diff_Rho(l,m,n),Diff_Phi(l,m,n)]=Cal_Difference_RhoPhi(P_rho,P_phi,Exp_rho,Exp_phi);
        end
    end
end
delete(kernum);
