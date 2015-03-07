function b=MPSS(CpRo,La,t,T1,T0,SS,BZ)
ae=La(4)/CpRo(4);

b=sum(SS,2).*La(4).*(T1-T0)./(pi*ae*t).^0.5;
b(BZ~=1)=0;

