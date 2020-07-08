% In this file, the iterator tolerance changed into 0.0001 (1e-4)
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

%%
% figure(1)
% plot(dis_theo,np_theo(:,1))
% hold on
% plot(dis_theo,np_theo(:,2))
%% Load the experimental data
data_para=csvread('../Data/para_power_100_600_itera_small_toler_1em4.csv',5,0);
data_ortho=csvread('../Data/ortho_power_100_600_itera_small_toler_1em4.csv',5,0);

dis_numer=data_para(:,1)*1e-9+radius;
p_numer=[data_para(:,4),data_ortho(:,4)];

p0=Fun_PowerFreeCOMSOL(lamda,1,1);
np_numer=p_numer/p0;

%% Plot and compare
figure(2)
plot(dis_theo-radius,np_theo(:,1),'r-')
hold on
plot(dis_numer-radius,np_numer(:,1),'r--')
hold on
plot(dis_theo-radius,np_theo(:,2),'b-')
hold on
plot(dis_numer-radius,np_numer(:,2),'b--')
grid on