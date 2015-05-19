function [W1,W6,Wo,Wg,W7,QQ,QQwo]=Well_MKTBO(P,Pw,Won,Uf,Cp,PR,Cpin,CMP,KWOG,v,n)

wmu=fusion(Cp(Won(:,1)),PR.mup);

TW=KWOG.w(v(Won(:,1)))./CMP.Bw(v(Won(:,1)),2)./wmu;
TO=KWOG.o(v(Won(:,1)))./CMP.Bo(v(Won(:,1)),2)/PR.mu(2);
TG=KWOG.g(v(Won(:,1)))./CMP.Bg(v(Won(:,1)),2)/PR.mu(3);
TP=KWOG.w(v(Won(:,1)))./wmu.*Cp(Won(:,1));

Tiw=(1-Cpin)./CMP.Bw(v(Won(:,1)),2)./wmu;
Tip=Cpin./CMP.Bw(v(Won(:,1)),2)./wmu;

W6 = Won(:,2).*(TW.*(Uf==1)+(Tiw+Tip).*(Uf==-1));
Wo = Won(:,2).*TO.*(Uf==1);
Wg = Won(:,2).*TG.*(Uf==1);
W1 = CMP.Cg(Won(:,1)).*Wg + (CMP.Rs(Won(:,1),2).*CMP.Cor(Won(:,1)) + CMP.Co(Won(:,1))).*Wo + CMP.Cw(Won(:,1)).*W6;

W7=-Won(:,2).*CMP.Cw(Won(:,1)).*(TP.*(Uf==1)+Tip.*(Uf==-1));%

dP = Pw - P(Won(:,1));
qwo(:,3) = Wg.*dP;
qwo(:,2) = Wo.*dP;
qwo(:,1) = W6.*dP;
QQ = CMP.Cg(Won(:,1)).*qwo(:,3) + (CMP.Rs(Won(:,1),2).*CMP.Cor(Won(:,1)) + CMP.Co(Won(:,1))).*qwo(:,2) + CMP.Cw(Won(:,1)).*qwo(:,1);
nw = size(Pw,1);
QQ =  sparse(Won(:,1),ones(1,nw),QQ,n,1);
QQwo(:,1) =  sparse(Won(:,1),ones(1,nw),qwo(:,1),n,1);
QQwo(:,2) =  sparse(Won(:,1),ones(1,nw),qwo(:,2),n,1);
%QQwo(:,3) =  sparse(Won(:,1),ones(1,nw),qwo(:,3),n,1);
