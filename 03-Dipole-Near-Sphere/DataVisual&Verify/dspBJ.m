function [dspBJ]=dspBJ(l,x)
l=l+0.5;
spBJ1=sqrt(pi/(x*2))*besselj(l,x);
spBJ2=sqrt(pi/(x*2))*besselj(l-1,x);
spBJ3=sqrt(pi/(x*2))*besselj(l+1,x);
dspBJ=0.5*spBJ1+x*0.5*(spBJ2-spBJ3);
end