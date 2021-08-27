function patUp=Cal_Pattern_1DDipole(p,DinLayer,objective)
% objective=1 means the pattern in objective
% objective=0 means the pattern in far field sphere surface
	UpEFar=Cal_Field_Dipole(p,DinLayer);
	patUp=(abs(UpEFar(:,:,1)).^2+abs(UpEFar(:,:,2)).^2+abs(UpEFar(:,:,3)).^2)./abs(cos(DinLayer.thetaUp).^(objective));
end