%% To cut the data from the experiment for laterl handling
load('./Data/BFP_tot_single_1.mat')

Nor_BFP=BFP_tot;
Nor_BFP(58,50)=(Nor_BFP(58,51)+Nor_BFP(58,49))/2;

Left=zeros(180,180);
Right=zeros(180,180);
Up=zeros(180,180);
Dn=zeros(180,180);

temp=linspace(1,180,180);


nLeft=13*ones(180);
nRight=167*ones(180);
nUp=167*ones(180);
nDn=13*ones(180);

num_kx=167-13+1;
num_ky=167-13+1;

%% 
figure(1)
pcolor(Nor_BFP)
hold on
plot(nLeft,temp,'color','white')
hold on
plot(nRight,temp,'color','white')
hold on
plot(temp,nUp,'color','white')
hold on
plot(temp,nDn,'color','white')
hold on
colormap jet
shading interp