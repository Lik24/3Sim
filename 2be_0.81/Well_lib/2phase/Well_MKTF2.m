function [W1,W6,Wo,W7]=Well_MKTF2(Wf,Won,Uf,Cp,mu,Cpin,kfw,kfo,CMP,CMPF)

TW=kfw(Won)./CMP.Bw(Won,2)/mu(1);
TO=kfo(Won)./CMP.Bo(Won,2)/mu(2);
TP=kfw(Won)/mu(4).*Cp(Won);

Tiw=(1-Cpin)/mu(1);%Sat_cal(ones(size(Won,1),1).*(1-Cpin),1,1,as,aw)/mu(1);
Tip=Cpin/mu(4);%Sat_cal(ones(size(Won,1),1).*Cpin,1,1,as,aw)/mu(4);

Wo = Wf.*TO.*(Uf==1);
W6 = Wf.*(TW.*(Uf==1)+(Tiw+Tip).*(Uf==-1));
W1 = CMPF.Cw(Won).*W6 + Wo;
W7=Wf.*(TP.*(Uf==1)+Tip.*(Uf==-1));

