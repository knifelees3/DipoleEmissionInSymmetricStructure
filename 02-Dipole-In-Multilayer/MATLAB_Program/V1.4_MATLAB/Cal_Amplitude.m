function [AUpS,AUpP,BUpS,BUpP]=Cal_Amplitude(RSUp, RPUp, RSDn, RPDn,kz0,POSD)
    AUpS=zeros(3,3);
    AUpP=zeros(3,3);
    BUpS=zeros(3,3);
    BUpP=zeros(3,3);
    
    
    AUpS(1:2, 1:2) = RSDn*(RSUp * exp(-1i * kz0*POSD) + exp(1i * kz0 * POSD))/(1 - RSUp * RSDn) +exp(-1i * kz0 * POSD);

    BUpS(1:2, 1:2) = RSUp * (RSDn * exp(1i * kz0 * POSD) +exp(-1i * kz0 * POSD))/(1 - RSUp * RSDn);
    
    AUpP(:, 1:2) = RPDn * (RPUp * exp(-1i * kz0 * POSD) - exp(1i * kz0 * POSD)) / (1 - RPUp * RPDn) +exp(-1i * kz0 * POSD);

    BUpP(:, 1:2) = -RPUp * (RPDn * exp(1i * kz0 * POSD) - exp(-1i *kz0 * POSD)) / (1 - RPUp * RPDn);
        
    AUpP(:, 3) = RPDn * (RPUp * exp(-1i * kz0 * POSD) + exp(1i * kz0 * POSD)) / (1 - RPUp * RPDn) + exp(-1i * kz0 * POSD);

    BUpP(:, 3) = RPUp * (RPDn * exp(1i * kz0 * POSD) + exp(-1i * kz0 * POSD)) / (1 - RPUp * RPDn);
    
end