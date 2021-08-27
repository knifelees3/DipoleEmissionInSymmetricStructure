function ShowStructure(DinLayer)
Eplist=DinLayer.Eplist;
dl=DinLayer.dl;
num_layer=DinLayer.num_layer;
num_dl=num_layer-1;
POSD=DinLayer.POSD;
thick=abs(dl(end)-dl(1));
dl_new=zeros(num_dl+2,1);
dl_new(2:end-1)=dl;
dl_new(1)=dl(1)-thick;
dl_new(end)=dl(end)+thick;

num_plot=10;
x=1:1:num_plot;
figure()
for l=2:num_dl+1
    plot(x,dl_new(l)*ones(num_plot,1)*1e9,'Linewidth',1.5);
    hold on
end
plot(x,dl_new(1)*ones(num_plot,1)*1e9,'--','Linewidth',1.5);
hold on
plot(x,dl_new(num_dl+2)*ones(num_plot,1)*1e9,'--','Linewidth',1.5);
hold on
plot(x(num_plot/2),POSD*1e9,'RO','Linewidth',1.5);
axis  off;
for l=1:num_layer
    text(x(2),(dl_new(l+1)+dl_new(l))/2*1e9,['Refractive Index = ',num2str(sqrt(Eplist(l)),3)],'fontname','times new roman','fontsize',12);
end
text(5.5,POSD*1e9,'Dipole Position','fontname','times new roman','fontsize',12);
title({'Schematic of the structure',['Wavelength=',num2str(DinLayer.WL0*1e9,3),'nm'],...
    ['Distance to lower interface:',num2str(DinLayer.dis*1e9,3),'nm']})
end