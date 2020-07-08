global k index_s index_air radius lamda k0 alpha alpham
index_s=sqrt(-50);
lamda=1e-6;
k0=2*pi/lamda;
k=k0;
radius=1e-7;
alpha=pi*2*radius/lamda;
alpham=alpha*index_s;
r=1.1e-7;
h=1e-8;
num=100;
A=zeros(num,2);
for m=1:num
	rate1=0;
    rate2=0;
	r=r+h;
	for l=1:80
    a1=(index_s^2*spBJ(l,alpham)*dspBJ(l,alpha)-spBJ(l,alpha)*dspBJ(l,alpham));
    a2=(index_s^2*spBJ(l,alpham)*dspBH(l,alpha)-spBH(l,alpha)*dspBJ(l,alpham));   
	b1=(spBJ(l,alpham)*dspBJ(l,alpha)-spBJ(l,alpha)*dspBJ(l,alpham));
    b2=(spBJ(l,alpham)*dspBH(l,alpha)-spBH(l,alpha)*dspBJ(l,alpham));
    x=k*r;
    al=a1/a2;
    bl=b1/b2;
	rate1=rate1+(2*l+1)*l*(l+1)*(-al)*(spBH(l,x)/(x))^2;
    rate2=rate2+(l+0.5)*((-al)*(dspBH(l,x)/(x))^2+(-bl)*(spBH(l,x))^2);   
end
rate1=real(rate1)*3/2+1;
A1(m,1)=r;
A1(m,2)=rate1;
rate2=1+1.5*real(rate2);
    A2(m,1)=r;
A2(m,2)=rate2;

 
end
A3=zeros(20,2);
A4=zeros(20,2);
A3(:,1)=per(:,1);
A3(:,2)=per(:,2)/0.042692;
A4(:,1)=para(:,1);
A4(:,2)=para(:,2)/0.042692;

figure(1)
plot(A1(:,1),A1(:,2),'s');
hold on;
plot(A3(:,1),A3(:,2),'-');
title('Perpendicular');
legend('theory','comsol')

figure(2);
plot(A2(:,1),A2(:,2),'s');
hold on;
plot(A4(:,1),A4(:,2),'-');
legend('theory','comsol');
title('Parallel');

