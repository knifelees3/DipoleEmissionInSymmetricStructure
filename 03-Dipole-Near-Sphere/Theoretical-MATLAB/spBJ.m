function [spBJ]=spBJ(l,x)
l=l+0.5;
spBJ=sqrt(pi/(x*2))*besselj(l,x);
end