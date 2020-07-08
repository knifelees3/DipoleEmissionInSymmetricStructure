index_s=sqrt(-22.473-1.3974i);
lamda=780e-9;
k0=2*pi/lamda;
k=k0;
radius=75e-9;
alpha=pi*2*radius/lamda;
alpham=alpha*index_s;
num_dis=100;
num_sum=80;
dis_theo=linspace(radius+100e-9,radius+600e-9,num_dis);
np_theo=Fun_nP_Cal(num_sum,num_dis,dis_theo,alpha,alpham,index_s,k);

%% Load the data from fdtd
data_ortho=csvread('../Data/purcell_ortho_fdtd.csv',0,0);
data_para=csvread('../Data/purcell_para_fdtd.csv',0,0);
p_fdtd_ortho=data_ortho(:,2);
p_fdtd_para=data_para(:,2);
dis_numer=data_ortho(:,1)+radius;
%% Load the data from meep
npt_ortho_mp=textread("../Data/npt_ortho_mp.txt");
dis=linspace(100e-9,600e-9,13);
%% Plot and compare
figure(2)
plot(dis_theo-radius,np_theo(:,1),'r-')
hold on
plot(dis_theo-radius,np_theo(:,2),'b-')
hold on
plot(dis_numer-radius,p_fdtd_para,'r*')
hold on
plot(dis_numer-radius,p_fdtd_ortho,'b*')
hold on
% plot(dis_numer-radius,p_fdtd_para,'r*')
% hold on
plot(dis,npt_ortho_mp(:,1),'b>')
legend('theo para','theo ortho','numer para lum','numer ortho lum','num ortho meep')