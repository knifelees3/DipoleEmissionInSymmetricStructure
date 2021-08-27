% Calculate the Green tensor for give single kx and ky

%% Calculate the RSP
function [GreenSUp,GreenPUp]=Cal_Green_Single(num_dl, Eplist, dl, klz, nUp, nDn, kx, ky,POSD,dUpFar)
num_layer=nUp+nDn+1;
[RSUp, RPUp, RSDn, RPDn, RS21, RP21, ~, ~] =Cal_RSP_kxy_Single(num_dl,Eplist, dl, klz, nUp, nDn);
%% Calculate the upper layers' green tensor
% The upper layer's green tensor
SSUpA=zeros(3,3);


PPUpA=zeros(3,3);

[PPUpA(:,:),~,SSUpA(:,:),~]=Cal_Green_1(kx,ky,klz(num_layer));


%% Calculate the upper amplitude of the emitting layer
AUpS0=zeros(3,3);
AUpP0=zeros(3,3);
BUpS0=zeros(3,3);
BUpP0=zeros(3,3);



[AUpS0(:,:),AUpP0(:,:),BUpS0(:,:),BUpP0(:,:)]=Cal_Amplitude(RSUp, RPUp, RSDn, RPDn,klz(nDn+1),POSD);

%% Cal culate the amplitude in the upper layer
AUpS=zeros(3,3,nUp+1);
AUpP=zeros(3,3,nUp+1);


[AUpS(:,:,:),AUpP(:,:,:),~,~]=Cal_Amplitude_Far(AUpS0(:,:),AUpP0(:,:),BUpS0(:,:),BUpP0(:,:),RS21(nDn+1:num_dl),RP21(nDn+1:num_dl),...
     klz(nDn+1:num_layer),dl(nDn+1:num_dl),nUp);     

 %% The Total Green tensor
 GreenSUp=zeros(3,3);
 GreenPUp=zeros(3,3);

 GreenSUp(:,:)=SSUpA(:,:).*AUpS(:,:,nUp+1).*exp(1i*klz(num_layer)*dUpFar);
 GreenPUp(:,:)=PPUpA(:,:).*AUpP(:,:,nUp+1).*exp(1i*klz(num_layer)*dUpFar);