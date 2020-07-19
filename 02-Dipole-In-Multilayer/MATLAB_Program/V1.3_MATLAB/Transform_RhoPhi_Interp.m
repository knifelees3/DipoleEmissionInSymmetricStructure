% Transform between different coordinate by using interplolation
% The input need the orginal gird and interpolated grid and the
% corresponding krho mat is needed for integration

function [Pattern_rho,Pattern_phi]=Transform_RhoPhi_Interp(kx_grid,ky_grid,kx_grid_in,ky_grid_in,PatternXY,krho)
Pattern_interp = interp2(kx_grid,ky_grid,PatternXY,kx_grid_in,ky_grid_in);
Pattern_rho=sum(Pattern_interp,1)/2/pi;
Pattern_phi=sum(Pattern_interp.*krho,2);
end