function Wf=KWell(K,H,S,Won)
Perf=0.5;
r0=0.05;
Sp=sum(S(Won,:),2);
Rk=(Sp/pi).^0.5;
Wcof=K(Won).*H(Won)./log(Rk/r0)*Perf;
Wf=2*pi*Wcof;