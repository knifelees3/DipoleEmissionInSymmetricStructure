function [E]=Cal_Field_Dipole(p,GreenSUp,GreenPUp,theta)
    [numkx, numky]=size(theta);
	ESUp=zeros(numkx,numky,3);
	EPUp=zeros(numkx,numky,3);
    
	for l=1:3
		ESUp(:,:,l)=GreenSUp(l,1,:,:)*p(1)+GreenSUp(l,2,:,:)*p(2)+GreenSUp(l,3,:,:)*p(3);
		EPUp(:,:,l)=GreenPUp(l,1,:,:)*p(1)+GreenPUp(l,2,:,:)*p(2)+GreenPUp(l,3,:,:)*p(3);
        ESUp(:,:,l)=ESUp(:,:,l).*(cos(theta));
        EPUp(:,:,l)=EPUp(:,:,l).*(cos(theta));
    end
    E=ESUp+EPUp;
end
