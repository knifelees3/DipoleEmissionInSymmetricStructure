function [AUpS,AUpP,BUpS,BUpP]=Cal_Amplitude_Far(AUpS0,AUpP0,BUpS0,BUpP0,RS21,RP21,klz,dl,nUp)

    AUpS=zeros(3,3,nUp+1);
    AUpP=zeros(3,3,nUp+1);
    BUpS=zeros(3,3,nUp+1);
    BUpP=zeros(3,3,nUp+1);
    
    AUpS(:,:,1)=AUpS0;
    AUpP(:,:,1)=AUpP0;
    BUpS(:,:,1)=BUpS0;
    BUpP(:,:,1)=BUpP0;
    

    for l=1:nUp

        expAA = exp(+1i * (klz(l) - klz(l + 1)) * dl(l));
        expAB = exp(-1i * (klz(l) + klz(l + 1)) * dl(l));
        expBA = exp(1i * (klz(l) + klz(l + 1)) * dl(l));
        expBB = exp(-1i * (klz(l) - klz(l + 1)) * dl(l));

        AUpS(:, :, l + 1) = 1 / (1 - RS21(l)) * (AUpS(:, :, l)* expAA + RS21(l) * BUpS(:, :, l) * expAB);
        BUpS(:, :, l + 1) = 1 / (1 - RS21(l)) * (RS21(l) * AUpS(:, :, l)* expBA + BUpS(:, :, l) * expBB);
        
        
        AUpP(:, 1:2, l + 1) = 1 / (1 + RP21(l)) * (AUpP(:, 1:2, l)* expAA + RP21(l) * BUpP(:, 1:2, l) * expAB); 
        BUpP(:, 1:2, l + 1) = 1 / (1 + RP21(l)) * (RP21(l)* AUpP(:, 1:2, l) * expBA + BUpP(:, 1:2, l) * expBB); 

        AUpP(:, 3, l + 1) = klz(l + 1) / klz(l) / (1 + RP21(l)) * (AUpP(:, 3, l) * expAA + RP21(l) * BUpP(:, 3, l) * expAB);    
        BUpP(:, 3, l + 1) = klz(l + 1) / klz(l) / (1 + RP21(l)) * (RP21(l) * AUpP(:, 3, l) * expBA + BUpP(:, 3, l) * expBB);    

%         AUpP(1:2, :, l + 1) = 1 / (1 + RP21(l)) * (AUpP(1:2, :, l)* expAA + RP21(l) * BUpP(1:2, :, l) * expAB); 
%         BUpP(1:2, :, l + 1) = 1 / (1 + RP21(l)) * (RP21(l)* AUpP(1:2, :, l) * expBA + BUpP(1:2, :, l) * expBB); 
% 
%         AUpP(3, :, l + 1) = klz(l + 1) / klz(l) / (1 + RP21(l)) * (AUpP(3, :, l) * expAA + RP21(l) * BUpP(3, :, l) * expAB);    
%         BUpP(3, :, l + 1) = klz(l + 1) / klz(l) / (1 + RP21(l)) * (RP21(l) * AUpP(3, :, l) * expBA + BUpP(3, :, l) * expBB);    
    end
end