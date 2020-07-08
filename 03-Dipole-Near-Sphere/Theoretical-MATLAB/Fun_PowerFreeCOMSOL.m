%% To calculate the emission power for a give dipole and material
% This is the function to be used compare with COMSOL
function P=Fun_PowerFreeCOMSOL(wavelength,epsilon,J)
	epsilon0 = 8.854187817620389850537e-12;
	mu0 = 1.2566370614359172953851e-6;
	c_const=1/sqrt(epsilon0*mu0);


	omega=2*pi./wavelength*c_const;
	n=sqrt(epsilon);

	% In COMSOL, the dipole moment is defined as current
	d=J*i/omega;

	P=(abs(d).^2).*omega.^4/3/c_const^3/4/pi/epsilon0/epsilon*n^3;
end