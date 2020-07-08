% To compare the iteration method when changing mesh and tolerance
% In this file, the iterator tolerance changed from 1e-3 to 1e-5
index_s = sqrt(-22.473 - 1.3974i);
lamda = 780e-9;
k0 = 2 * pi / lamda;
k = k0;
radius = 75e-9;
alpha = pi * 2 * radius / lamda;
alpham = alpha * index_s;
num_dis = 100;
num_sum = 80;
dis_theo = linspace(radius + 100e-9, radius + 600e-9, num_dis);
np_theo = Fun_nP_Cal(num_sum, num_dis, dis_theo, alpha, alpham, index_s, k);

%%
% figure(1)
% plot(dis_theo,np_theo(:,1))
% hold on
% plot(dis_theo,np_theo(:,2))
%% Load the experimental data
data_para_1 = csvread('../Data/para_power_100_600_itera_small_toler_1em3.csv', 5, 0);
data_ortho_1 = csvread('../Data/ortho_power_100_600_itera_small_toler_1em3.csv', 5, 0);
data_para_2 = csvread('../Data/para_power_100_600_itera_small_toler_1em4.csv', 5, 0);
data_ortho_2 = csvread('../Data/ortho_power_100_600_itera_small_toler_1em4.csv', 5, 0);
data_para_3 = csvread('../Data/para_power_100_600_itera_small_toler_1em5.csv', 5, 0);
data_ortho_3 = csvread('../Data/ortho_power_100_600_itera_small_toler_1em5.csv', 5, 0);

dis_numer = data_para(:, 1) * 1e-9 + radius;
p_numer = [data_para_1(:, 4), data_ortho_1(:, 4), data_para_2(:, 4), data_ortho_2(:, 4), data_para_3(:, 4), data_ortho_3(:, 4)];

p0 = Fun_PowerFreeCOMSOL(lamda, 1, 1);
np_numer = p_numer / p0;

%% Plot and compare

dis_numer_re = (dis_numer - radius) * 1e9;
dis_theo_re = (dis_theo - radius) * 1e9;
figure(2)
plot(dis_theo_re, np_theo(:, 1), 'r--','linewidth',1)
hold on
plot(dis_theo_re, np_theo(:, 2), 'b--','linewidth',1)
hold on
plot(dis_numer_re, np_numer(:, 1), '-');%'markersize',10)
hold on
plot(dis_numer_re, np_numer(:, 2), '-');%'markersize',10)
hold on
plot(dis_numer_re, np_numer(:, 3), '-');%'markersize',10)
hold on
plot(dis_numer_re, np_numer(:, 4), '-');%'markersize',10)
hold on
plot(dis_numer_re, np_numer(:, 5), '-');%'markersize',10)
hold on
plot(dis_numer_re, np_numer(:, 6), '-');%'markersize',10)
set(gca,'fontname','times new roman', 'fontsize',15,'linewidth',1.5)
xlim([100, 600]);
xlabel('Distance (nm)')
ylabel('\Gamma_{t}/\Gamma_{0}')

%% Plot with error bar
