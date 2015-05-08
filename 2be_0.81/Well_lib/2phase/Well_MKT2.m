function [W1,W6,Wo,W7]=Well_MKT2(Won,Uf,Cp,PR,Cpin,CMP,KWOG,v)

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

W6 = Won(:,2).*(TW.*(Uf==1)+(Tiw+Tip).*(Uf==-1));
Wo = Won(:,2).*TO.*(Uf==1);
W7=-Won(:,2).*(TP.*(Uf==1)+Tip.*(Uf==-1));%
W1 = CMP.Cw(v(Won(:,1))).*W6 + Wo;


