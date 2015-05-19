function [W1,W6,Wo,Wg,W6S,WoSo,WoSw,WgS,W7]=Well_MKTSeqBO(Won,Uf,Cp,PR,Cpin,CMP,KWOG,v,Pw,P)

wmu=fusion(Cp(Won(:,1)),PR.mup);
dP = Pw(Won(:,3)) - P(Won(:,1));
TWS = dP./CMP.Bw(v(Won(:,1)),2)./wmu.*KWOG.dwds(v(Won(:,1)));
TOSO = -dP.*KWOG.dodsg(v(Won(:,1)))./CMP.Bo(v(Won(:,1)),2)/PR.mu(2);
TOSW = dP.*(KWOG.dodsw(v(Won(:,1)))-KWOG.dodsg(v(Won(:,1))))./CMP.Bo(v(Won(:,1)),2)/PR.mu(2);
TGS = -dP.*KWOG.dgds(v(Won(:,1)))./CMP.Bo(v(Won(:,1)),2)/PR.mu(3);

TW=KWOG.w(v(Won(:,1)))./CMP.Bw(v(Won(:,1)),2)./wmu;
TO=KWOG.o(v(Won(:,1)))./CMP.Bo(v(Won(:,1)),2)/PR.mu(2);
TG=KWOG.g(v(Won(:,1)))./CMP.Bg(v(Won(:,1)),2)/PR.mu(3);
TP=KWOG.w(v(Won(:,1)))./wmu.*Cp(Won(:,1));

Tiw=(1-Cpin)./CMP.Bw(v(Won(:,1)),2)./wmu;
Tip=Cpin./CMP.Bw(v(Won(:,1)),2)./wmu;

W6 = Won(:,2).*(TW.*(Uf==1)+(Tiw+Tip).*(Uf==-1));
W6S = Won(:,2).*TWS.*(Uf==1);
WgS = Won(:,2).*TGS.*(Uf==1);
Wo = Won(:,2).*TO.*(Uf==1);
WoSo = Won(:,2).*TOSO.*(Uf==1);
WoSw = Won(:,2).*TOSW.*(Uf==1);
Wg = Won(:,2).*TG.*(Uf==1);
W1 = CMP.Cg(v(Won(:,1))).*Wg + (CMP.Rs(v(Won(:,1)),2).*CMP.Cor(v(Won(:,1))) + CMP.Co(v(Won(:,1)))).*Wo + CMP.Cw(v(Won(:,1))).*W6;

W7=-Won(:,2).*CMP.Cw(v(Won(:,1))).*(TP.*(Uf==1)+Tip.*(Uf==-1));