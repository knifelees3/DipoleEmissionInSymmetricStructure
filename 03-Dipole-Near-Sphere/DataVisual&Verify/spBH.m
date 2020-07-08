function [spBH]=spBH(l,x)
l=l+0.5;
spBH=sqrt(pi/(x*2))*(besselj(l,x)+1j*bessely(l,x));
end