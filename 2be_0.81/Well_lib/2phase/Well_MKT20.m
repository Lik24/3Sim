function [W1,W6,Wo,W7,QQ,QQw]=Well_MKT20(P,Pw,Won,Uf,Cp,aw,as,PR,Cpin,CMP,KWOG,v,n)

wmu=fusion(Cp(Won(:,1)),PR.mup);

TW=KWOG.w(v(Won(:,1)))./CMP.Bw(v(Won(:,1)),2)./wmu;
TO=KWOG.o(v(Won(:,1)))./CMP.Bo(v(Won(:,1)),2)/PR.mu(2);
TP=KWOG.w(v(Won(:,1)))./wmu.*Cp(Won(:,1));
% 
% TW=kfw(Won)/mu(1);
% TO=kfo(Won),2,1,as,aw)/mu(2);
% TP=Sat_cal(Sw(Won),1,1,as,aw)/mu(4).*Cp(Won);

Tiw=(1-Cpin)./wmu;
Tip=Cpin./wmu;

% size((TW+TO+TP))
% size((Uf==1))
% size((Tiw+Tip))
% size((Uf==-1))
% size(Wf)
% size(TW)
% size(Uf)

W6 = Won(:,2).*(TW.*(Uf==1)+(Tiw+Tip).*(Uf==-1));
Wo = Won(:,2).*TO.*(Uf==1);
W1 = CMP.Cw(v(Won(:,1))).*W6 + Wo;

W7=-Won(:,2).*(TP.*(Uf==1)+Tip.*(Uf==-1));%

dP = Pw - P(Won(:,1));
qwo(:,2) = Wo.*dP;
qwo(:,1) = W6.*dP;
QQ = CMP.Cw(v(Won(:,1))).*qwo(:,1) + qwo(:,2);
nw = size(Pw,1);
QQ =  sparse(Won(:,1),ones(1,nw),QQ,n,1);
QQw(:,1) =  sparse(Won(:,1),ones(1,nw),qwo(:,1),n,1);
QQw(:,2) =  sparse(Won(:,1),ones(1,nw),qwo(:,2),n,1);

 

