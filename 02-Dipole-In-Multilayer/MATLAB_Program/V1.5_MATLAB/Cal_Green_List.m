
function DinLayer=Cal_Green_List(DinLayer)
% This program returns the value of single kx and ky

% Extract some paremeters
num_kx=DinLayer.num_kx;
num_ky=DinLayer.num_ky;
nUp=DinLayer.nUp;
nDn=DinLayer.nDn;
dUpFar=DinLayer.dUpFar;
dDnFar=DinLayer.dDnFar;
kx_grid=DinLayer.kx_grid;
ky_grid=DinLayer.ky_grid;
num_dl=DinLayer.num_dl;
dl=DinLayer.dl;
klz=DinLayer.klz;
Eplist=DinLayer.Eplist;
num_layer=DinLayer.num_layer;
POSD=DinLayer.POSD;
%% Calculate the RSP
[RSUp, RPUp, RSDn, RPDn, RS21, RP21, RS12, RP12] =Cal_RSP_kxy(num_dl, num_kx, num_ky, Eplist, dl, klz, nUp, nDn);

%% Calculate the upper layers' green tensor
% The upper layer's green tensor
SSUpA=zeros(3,3,num_kx,num_ky);
SSUpB=zeros(3,3,num_kx,num_ky);

PPUpA=zeros(3,3,num_kx,num_ky);
PPUpB=zeros(3,3,num_kx,num_ky);

for l=1:num_kx
    for m=1:num_ky
        [PPUpA(:,:,l,m),PPUpB(:,:,l,m),SSUpA(:,:,l,m),SSUpB(:,:,l,m)]=Cal_Green_Step1(kx_grid(l,m),ky_grid(l,m),DinLayer.klz(l,m,num_layer));
    end
end

%% Calculate the upper amplitude of the emitting layer
AUpS0=zeros(3,3,num_kx,num_ky);
AUpP0=zeros(3,3,num_kx,num_ky);
BUpS0=zeros(3,3,num_kx,num_ky);
BUpP0=zeros(3,3,num_kx,num_ky);


for l=1:num_kx
    for m=1:num_ky
    [AUpS0(:,:,l,m),AUpP0(:,:,l,m),BUpS0(:,:,l,m),BUpP0(:,:,l,m)]=Cal_Amplitude(RSUp(l,m), RPUp(l,m), RSDn(l,m), RPDn(l,m),klz(l,m,nDn+1),POSD);
    end
end

%% Cal culate the amplitude in the upper layer
AUpS=zeros(3,3,num_kx,num_ky,nUp+1);
AUpP=zeros(3,3,num_kx,num_ky,nUp+1);
BUpS=zeros(3,3,num_kx,num_ky,nUp+1);
BUpP=zeros(3,3,num_kx,num_ky,nUp+1);
 for l=1:num_kx
     for m=1:num_ky
         [AUpS(:,:,l,m,:),AUpP(:,:,l,m,:),BUpS(:,:,l,m,:),BUpP(:,:,l,m,:)]=...
             Cal_Amplitude_Far(AUpS0(:,:, l, m),AUpP0(:,:,l,m),BUpS0(:,:,l,m),BUpP0(:,:,l,m),...
             RS21(l,m,nDn+1:num_dl),RP21(l,m,nDn+1:num_dl),...
             klz(l,m,nDn+1:num_layer),dl(nDn+1:num_dl),nUp);     
     end
 end
  
 %% The Total Green tensor
 GreenSUp=zeros(3,3,num_kx,num_ky);
 GreenPUp=zeros(3,3,num_kx,num_ky);
 for l=1:num_kx
     for m=1:num_ky
 GreenSUp(:,:,l,m)=SSUpA(:,:,l,m).*AUpS(:,:,l,m,nUp+1).*exp(1i*klz(l,m,num_layer)*dUpFar);
 GreenPUp(:,:,l,m)=PPUpA(:,:,l,m).*AUpP(:,:,l,m,nUp+1).*exp(1i*klz(l,m,num_layer)*dUpFar);
     end
 end
 
DinLayer.GreenSUp=GreenSUp;
DinLayer.GreenPUp=GreenPUp;
% *******************************************
showtext=strcat(datestr(now,	'yyyy-mm-dd HH:MM:SS'),': Green Tensor Is Calculated \n');
fprintf(showtext);
end