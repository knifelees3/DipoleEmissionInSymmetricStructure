function DinLayerRever=StructureReSet(DinLayer)
	DinLayerRever=DinLayer;
	DinLayerRever.nUp=DinLayer.nDn;
	DinLayerRever.nDn=DinLayer.nUp;

    
	for l=1:DinLayerRever.num_dl
	DinLayerRever.dl(DinLayerRever.num_dl-l+1) = -DinLayer.dl(l);
    end
    for l=1:DinLayerRever.num_layer
    DinLayerRever.Eplist(DinLayerRever.num_layer-l+1) = DinLayer.Eplist(l);
    end
    
	DinLayerRever.POSD = DinLayerRever.dl(DinLayerRever.nDn+1) - DinLayerRever.dis;

    DinLayerRever.kl = DinLayerRever.k0 * sqrt(DinLayerRever.Eplist);
    DinLayerRever.ke=DinLayerRever.k0*sqrt(DinLayerRever.Eplist(DinLayerRever.num_layer));
    

	% The Z component of the wave vector
    DinLayerRever.klz = zeros(DinLayerRever.num_kx, DinLayerRever.num_ky, DinLayerRever.num_layer);
    DinLayerRever.theta = zeros(DinLayerRever.num_kx, DinLayerRever.num_ky, DinLayerRever.num_layer);
    DinLayerRever.krho_grid=sqrt(DinLayerRever.kx_grid.^2+DinLayerRever.ky_grid.^2);
    for l =1:DinLayerRever.num_layer
        DinLayerRever.klz(:, :, l) = sqrt(DinLayerRever.kl(l).^2 - DinLayerRever.krho_grid.^2);
        DinLayerRever.theta(:, :, l) = asin(DinLayerRever.krho_grid./DinLayerRever.kl(l));
    end

	DinLayerRever.thetaUp=DinLayerRever.theta(:,:,DinLayerRever.num_layer);
end