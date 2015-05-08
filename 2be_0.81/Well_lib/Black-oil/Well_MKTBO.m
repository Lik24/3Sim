function [W1,W6,Wo,Wg,W7]=Well_MKT(Won,Uf,Cp,PR,Cpin,CMP,KWOG,v)

wmu=fusion(Cp(Won(:,1)),PR.mup);

TW=KWOG.w(v(Won(:,1)))./CMP.Bw(v(Won(:,1)),2)./wmu;
TO=KWOG.o(v(Won(:,1)))./CMP.Bo(v(Won(:,1)),2)/PR.mu(2);
TG=KWOG.g(v(Won(:,1)))./CMP.Bg(v(Won(:,1)),2)/PR.mu(3);
TP=KWOG.w(v(Won(:,1)))./wmu.*Cp(Won(:,1));

Tiw=(1-Cpin)./wmu;
Tip=Cpin./wmu;

W6 = Won(:,2).*(TW.*(Uf==1)+(Tiw+Tip).*(Uf==-1));
Wo = Won(:,2).*TO.*(Uf==1);
Wg = Won(:,2).*TG.*(Uf==1);
W1 = CMP.Cg(v(Won(:,1))).*Wg + (CMP.Rs(v(Won(:,1)),2).*CMP.Cor(v(Won(:,1))) + CMP.Co(v(Won(:,1)))).*Wo + CMP.Cw(v(Won(:,1))).*W6;

W7=-Won(:,2).*CMP.Cw(v(Won(:,1))).*(TP.*(Uf==1)+Tip.*(Uf==-1));

