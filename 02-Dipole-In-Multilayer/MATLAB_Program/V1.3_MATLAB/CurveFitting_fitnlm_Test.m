rand_data =random('norm',2,0.3,100,1);

Gau_Fun=@(b,x)b(1)+b(2)*exp(-(x-b(3)).^2);
beta0=[1,1,1];

x=linspace(0,10,100);
ba=[2,2,3];
Exp_1=Gau_Fun(ba,x)+rand_data'/2;

beta_guess=[1,2,3.3];
mdl=fitnlm(x',Exp_1,Gau_Fun,beta_guess);
coe=mdl.Coefficients.Estimate;
Exp_2=Gau_Fun(coe,x);
figure(1)
plot(x,Exp_1,'*')
%hold on
%plot(x,Exp_1-rand_data'/3,'*')
hold on
plot(x,Exp_2,'-')
legend('Original','Fitted')


