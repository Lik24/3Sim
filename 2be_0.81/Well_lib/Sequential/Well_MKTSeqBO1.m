function [W1,W6,Wo,Wg,W6S,WoSo,WoSw,W7,qq,qqwog]=Well_MKTSeqBO1(P,Pw,Won,Uf,Cp,PR,Cpin,CMP,KWOG,v,n)

wmu=fusion(Cp(Won(:,1)),PR.mup);

TiwS = (Pw(Won(:,3))-P(Won(:,1)))./CMP.Bw(v(Won(:,1)),2)./wmu;
TWS = TiwS.*KWOG.dwds(v(Won(:,1)));
TOSO = -(Pw(Won(:,3))-P(Won(:,1))).*KWOG.dodsg(v(Won(:,1)))./CMP.Bo(v(Won(:,1)),2)/PR.mu(2);
TOSW = (Pw(Won(:,3))-P(Won(:,1))).*(KWOG.dodsw(v(Won(:,1)))-KWOG.dodsg(v(Won(:,1))))./CMP.Bo(v(Won(:,1)),2)/PR.mu(2);

TW=KWOG.w(v(Won(:,1)))./CMP.Bw(v(Won(:,1)),2)./wmu;
TO=KWOG.o(v(Won(:,1)))./CMP.Bo(v(Won(:,1)),2)/PR.mu(2);
TG=KWOG.g(v(Won(:,1)))./CMP.Bg(v(Won(:,1)),2)/PR.mu(3);
TP=KWOG.w(v(Won(:,1)))./wmu.*Cp(Won(:,1));

Tiw=(1-Cpin)./CMP.Bw(v(Won(:,1)),2)./wmu;
Tip=Cpin./CMP.Bw(v(Won(:,1)),2)./wmu;

W6 = Won(:,2).*(TW.*(Uf==1)+(Tiw+Tip).*(Uf==-1));
W6S = Won(:,2).*(TWS.*(Uf==1)+TiwS.*(Uf==-1));
Wo = Won(:,2).*TO.*(Uf==1);
WoSo = Won(:,2).*TOSO.*(Uf==1);
WoSw = Won(:,2).*TOSW.*(Uf==1);
Wg = Won(:,2).*TG.*(Uf==1);
W1 = CMP.Cg(v(Won(:,1))).*Wg + (CMP.Rs(v(Won(:,1)),2).*CMP.Cor(v(Won(:,1))) + CMP.Co(v(Won(:,1)))).*Wo + CMP.Cw(v(Won(:,1))).*W6;

W7=-Won(:,2).*CMP.Cw(v(Won(:,1))).*(TP.*(Uf==1)+Tip.*(Uf==-1));

dP = Pw - P(Won(:,1));
qwog(:,3) = Wg.*dP;
qwog(:,2) = Wo.*dP;
qwog(:,1) = W6.*dP;
qq = W1.*dP;
nw = size(Pw,1);
qq =  sparse(Won(:,1),ones(1,nw),qq,n,1);
qqwog(:,1) =  sparse(Won(:,1),ones(1,nw),qwog(:,1),n,1);
qqwog(:,2) =  sparse(Won(:,1),ones(1,nw),qwog(:,2),n,1);
qqwog(:,3) =  sparse(Won(:,1),ones(1,nw),qwog(:,3),n,1);

