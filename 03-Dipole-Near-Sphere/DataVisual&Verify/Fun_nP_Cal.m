function nP=Fun_nP_Cal(num_sum,num_dis,dis,alpha,alpham,index_s,k)
nP=zeros(num_dis,2); % with nP(:,1) store the np_ortho and nP(:,2) store the np_para
for m=1:num_dis
	rate_ortho=0;
    rate_para=0;
    r=dis(m);
	for l=1:num_sum
        a1=(index_s^2*spBJ(l,alpham)*dspBJ(l,alpha)-spBJ(l,alpha)*dspBJ(l,alpham));
        a2=(index_s^2*spBJ(l,alpham)*dspBH(l,alpha)-spBH(l,alpha)*dspBJ(l,alpham));   
        b1=(spBJ(l,alpham)*dspBJ(l,alpha)-spBJ(l,alpha)*dspBJ(l,alpham));
        b2=(spBJ(l,alpham)*dspBH(l,alpha)-spBH(l,alpha)*dspBJ(l,alpham));
	    x=k*r;
	    al=a1/a2;
	    bl=b1/b2;
		rate_ortho=rate_ortho+(2*l+1)*l*(l+1)*(-al)*(spBH(l,x)/(x))^2;
	    rate_para=rate_para+(l+0.5)*((-al)*(dspBH(l,x)/(x))^2+(-bl)*(spBH(l,x))^2);   
	end
	rate_ortho=real(rate_ortho)*3/2+1;
	rate_para=1+1.5*real(rate_para);
	nP(m,1)=rate_ortho;
	nP(m,2)=rate_para; 
end
end