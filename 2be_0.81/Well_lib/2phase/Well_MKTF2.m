function [W1,W6,Wo,W7]=Well_MKTF2(Wf,Won,Uf,Cp,mu,Cpin,kfw,kfo,CMP,v)

TW=kfw(Won)./CMP.Bw(v(Won),2)/mu(1);
TO=kfo(Won)./CMP.Bo(v(Won),2)/mu(2);
TP=kfw(Won)/mu(4).*Cp(Won);

Tiw=(1-Cpin)/mu(1);
Tip=Cpin/mu(4);

Wo = Wf.*TO.*(Uf==1);
W6 = Wf.*(TW.*(Uf==1)+(Tiw+Tip).*(Uf==-1));
W1 = CMP.Cw(v(Won)).*W6 + Wo;
W7=Wf.*(TP.*(Uf==1)+Tip.*(Uf==-1));

