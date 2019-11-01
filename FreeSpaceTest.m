%% Load the exparimental data
Numerical=csvread('pt_freespace_700_900.csv',5,0);

% Theoretical part

wavelength=Numerical(:,1);
epsilon0 = 8.854187817620389850537e-12;
mu0 = 1.2566370614359172953851e-6;
c_const=1/sqrt(epsilon0*mu0);
omega=2*pi./wavelength*c_const;
% in COMSOL, the dipole moment is defined with current density
d=1i./omega;
P_theo=power_free(wavelength,1,d);

P_simu=Numerical(:,3:5);
%%
figure(1)
plot(wavelength,P_theo,'-')
hold on
plot(wavelength,P_simu(:,1),'s')
hold on
plot(wavelength,P_simu(:,2),'s')
hold on
plot(wavelength,P_simu(:,3),'s')