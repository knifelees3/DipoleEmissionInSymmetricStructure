function Pattern=Cal_Pattern_QD_single(p1,p2,GreenSUp,GreenPUp,theta)
	ESUp1=zeros(3,1);
	EPUp1=zeros(3,1);
	ESUp2=zeros(3,1);
	EPUp2=zeros(3,1);

	for l=1:3
		ESUp1(l)=GreenSUp(l,1)*p1(1)+GreenSUp(l,2)*p1(2)+GreenSUp(l,3)*p1(3);
		EPUp1(l)=GreenPUp(l,1)*p1(1)+GreenPUp(l,2)*p1(2)+GreenPUp(l,3)*p1(3);
		ESUp2(l)=GreenSUp(l,1)*p2(1)+GreenSUp(l,2)*p2(2)+GreenSUp(l,3)*p2(3);
		EPUp2(l)=GreenPUp(l,1)*p2(1)+GreenPUp(l,2)*p2(2)+GreenPUp(l,3)*p2(3);
	end
	% Cal the emission pattern 
	PatternS=(abs(ESUp1(1)).^2+abs(ESUp1(2)).^2+abs(ESUp1(3)).^2+...
	    abs(ESUp2(1)).^2+abs(ESUp2(2)).^2+abs(ESUp2(3)).^2).*abs(cos(theta)).^2;
	PatternP=(abs(EPUp1(1)).^2+abs(EPUp1(2)).^2+abs(EPUp1(3)).^2+...
        abs(EPUp2(1)).^2+abs(EPUp2(2)).^2+abs(EPUp2(3)).^2).*abs(cos(theta)).^2;
	Pattern=PatternS+PatternP;
end
