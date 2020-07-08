function [dspBH]=dspBH(l,x)
l=l+0.5;
spBH1=sqrt(pi/(x*2))*(besselj(l,x)+1j*bessely(l,x));
spBH2=sqrt(pi/(x*2))*(besselj(l-1,x)+1j*bessely(l-1,x));
spBH3=sqrt(pi/(x*2))*(besselj(l+1,x)+1j*bessely(l+1,x));
dspBH=0.5*spBH1+x*0.5*(spBH2-spBH3);
end