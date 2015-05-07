function [W1,W6,Wo,W7,QQ,QQw]=Well_MKTF20(P,Pw,Wf,Won,Uf,Cp,mu,Cpin,kfw,kfo,SGM,CMP)

TW=kfw(Won)./CMP.Bw(Won,2)/mu(1);
TO=kfo(Won)./CMP.Bo(Won,2)/mu(2);
TP=kfw(Won)/mu(4).*Cp(Won);

Tiw=(1-Cpin)/mu(1);%Sat_cal(ones(size(Won,1),1).*(1-Cpin),1,1,as,aw)/mu(1);
Tip=Cpin/mu(4);%Sat_cal(ones(size(Won,1),1).*Cpin,1,1,as,aw)/mu(4);

Wo = Wf.*TO.*(Uf==1);
W6 = Wf.*(TW.*(Uf==1)+(Tiw+Tip).*(Uf==-1));
W1 = CMP.Cw(Won).*W6 + Wo;
W7=Wf.*(TP.*(Uf==1)+Tip.*(Uf==-1));%

dP = Pw - P(Won(:,1));
qwo(:,2) = Wo.*dP;
qwo(:,1) = W6.*dP;
QQ = CMP.Cw(Won(:,1)).*qwo(:,1) + qwo(:,2);
nw = size(Pw,1);
n = size(P,1);
QQ =  sparse(Won(:,1),ones(1,nw),QQ,n,1);
QQw(:,1) =  sparse(Won(:,1),ones(1,nw),qwo(:,1),n,1);
QQw(:,2) =  sparse(Won(:,1),ones(1,nw),qwo(:,2),n,1);
