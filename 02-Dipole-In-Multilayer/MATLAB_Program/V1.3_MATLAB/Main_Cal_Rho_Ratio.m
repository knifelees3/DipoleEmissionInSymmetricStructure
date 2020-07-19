%% To calculate the distribution along rho and sweep the ratio of the BFP edge point and the maximum point
num_alpha=20;
num_phi_1=20;
num_phi_2=20;

alpha_mat=linspace(0,pi,num_alpha);
phi_1_mat=linspace(0,pi*2,num_phi_1);
phi_2_mat=linspace(0,pi*2,num_phi_2);


Ratio_rho=zeros(num_alpha,num_phi_1,num_phi_2);
Prho=zeros(num_alpha,num_phi_1,num_phi_2,num_kx);

% kernum=parpool(60);
for l=1:num_alpha
    %showtext=strcat(datestr(now,'yyyy-mm-dd HH:MM:SS'),': current loop is ',num2str(l),'\n');
    %fprintf(showtext);
    for m=1:num_phi_1
        for n=1:num_phi_2
            [p1,p2]=DipoleQD(alpha_mat(l),phi_1_mat(m),phi_2_mat(n));
            PatternMat=Cal_Pattern_QD(num_kx,num_ky,p1,p2,GreenSUp,GreenPUp,thetaUp);
            Prho(l,m,n,:)=sum(PatternMat,1)/max(sum(PatternMat,1));
            Ratio_rho(l,m,n)=Prho(l,m,n,num_kx);
        end
    end
end

% delete(kernum);

%%
% p1=[1,0,0]';
% p2=[1,0,0]';
% %[p1,p2]=DipoleQD(alpha_mat(l),phi_1_mat(m),phi_2_mat(n));
% PatternMat=Cal_Pattern_QD(num_kx,num_ky,p1,p2,GreenSUp,GreenPUp,thetaUp);
% 
% 
% figure(1)
% pcolor(kx_grid,ky_grid,PatternMat);shading interp;colormap jet;
% 
% P_rho=sum(PatternMat,1);
% P_phi=sum(PatternMat,2);
% 
% plot(kphi,P_phi)
%% Two special case
p1=[0,0,1]';
p2=[0,0,1]';
PatternMat=Cal_Pattern_QD(num_kx,num_ky,p1,p2,GreenSUp,GreenPUp,thetaUp);
Prho_1=sum(PatternMat,1)/max(sum(PatternMat,1));

p1=[1,0,0]';
p2=[1,0,0]';
PatternMat=Cal_Pattern_QD(num_kx,num_ky,p1,p2,GreenSUp,GreenPUp,thetaUp);
Prho_2=sum(PatternMat,1)/max(sum(PatternMat,1));

%%
figure(2)
for l=1:num_alpha
    for m=1:num_phi_1
        for n=1:num_phi_2
         plot(krho/k0,reshape(Prho(l,m,n,:),num_kx,1))
         hold on
        end
    end
end
plot(krho/k0,reshape(Prho_1,num_kx,1),'r*')
hold on
plot(krho/k0,reshape(Prho_2,num_kx,1),'b*')
hold on
xlabel('k_{\rho}/k0')