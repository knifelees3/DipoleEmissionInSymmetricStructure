p1=[0,0,1]';
p2=[0,0,0]';

Main_Basic_1;
Main_Green;

Pattern1=Cal_Pattern_QD(num_kx,num_ky,p1,p2,GreenSUp,GreenPUp,thetaUp);
nPattern1=Pattern1/max(max(Pattern1));
P_rho_1=sum(nPattern1,1).*sin(Theta_mat);

%%
Main_Basic_2;
Main_Green;

Pattern2=Cal_Pattern_QD(num_kx,num_ky,p1,p2,GreenSUp,GreenPUp,thetaUp);
nPattern2=Pattern2/max(max(Pattern1));
P_rho_2=sum(nPattern2,1).*sin(Theta_mat);


Main_Basic_3;
Main_Green;

Pattern3=Cal_Pattern_QD(num_kx,num_ky,p1,p2,GreenSUp,GreenPUp,thetaUp);
nPattern3=Pattern3/max(max(Pattern1));
P_rho_3=sum(nPattern3,1).*sin(Theta_mat);

%% Plot together
figure(1)
plot((Theta_mat)/pi*180,P_rho_1/max(P_rho_1),'g')
hold on
plot((Theta_mat)/pi*180,P_rho_2/max(P_rho_1),'r')
hold on
plot((Theta_mat)/pi*180,P_rho_3/max(P_rho_1)/1.78,'k')
hold on
grid on
xlim([0,90])
