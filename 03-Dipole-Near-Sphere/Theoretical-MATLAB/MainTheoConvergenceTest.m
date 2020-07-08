% To test when will the theoretical expressions will convengent
index_s=sqrt(-22.473-1.3974i);
%index_s=sqrt(-50);
lamda=1000e-9;
k0=2*pi/lamda;
k=k0;
radius=110e-9;
alpha=pi*2*radius/lamda;
alpham=alpha*index_s;
num_dis=100;

num_sum=[1,3,50];
num_len=length(num_sum);
dis_theo=linspace(radius+50e-9,radius+1000e-9,num_dis);

np_theo=zeros(num_dis,2,3);

for l=1:num_len
    np_theo(:,:,l)=Fun_nP_Cal(num_sum(l),num_dis,dis_theo,alpha,alpham,index_s,k);
end

%% 
figure(1)
for l=1:num_len
    plot(dis_theo,np_theo(:,1,l),'-')
    hold on
    plot(dis_theo,np_theo(:,2,l),'--')
end
