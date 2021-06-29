% To obatin the dipole far field pattern distribution

function  [pat,pat_rho,pat_phi]=Get_Pattern_3D_Dipole_Dark_Axis(paramat,numcut)
    numin=156;
    NARange=1.4/1.4579;
    uxFDTD=linspace(-NARange,NARange,numin);
    uyFDTD=linspace(-NARange,NARange,numin);
    
    [uxFDTDgrid,uyFDTDgrid]=meshgrid(uxFDTD,uyFDTD);
    urhoFDTDgrid=sqrt(uxFDTDgrid.^2+uyFDTDgrid.^2);
    thetaFDTD=asin(urhoFDTDgrid);

    thetamatin=linspace(0,pi/2,numin);
    phimatin=linspace(0,2*pi,numin);
    [thetagridin,phigridin]=meshgrid(thetamatin,phimatin);

    ux_grid_in=NARange*sin(thetagridin).*cos(phigridin);
    uy_grid_in=NARange*sin(thetagridin).*sin(phigridin);
    urho_grid_in=sqrt(ux_grid_in.^2+uy_grid_in.^2);
       
    [d1,d2,d3]=DipoleQDDark(paramat);
    
    if nargin > 1
            [pat1,~,~]=Get_Pattern_Coherent_Combination(d1,numcut);
            [pat2,~,~]=Get_Pattern_Coherent_Combination(d2,numcut);
            [pat3,~,~]=Get_Pattern_Coherent_Combination(d3,numcut);
    else
            [pat1,~,~]=Get_Pattern_Coherent_Combination(d1);
            [pat2,~,~]=Get_Pattern_Coherent_Combination(d2);
            [pat3,~,~]=Get_Pattern_Coherent_Combination(d3);
    end
    pat=pat1+pat2+pat3;
    [pat_rho,pat_phi]=Transform_RhoPhi_Interp(...
    uxFDTD,uyFDTD,ux_grid_in,uy_grid_in,pat,urho_grid_in);
end