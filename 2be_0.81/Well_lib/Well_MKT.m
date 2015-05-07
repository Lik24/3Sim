function [W1,W6,Wo,Wg,W7,QQ]=Well_MKT(dP,dPw,Won,Uf,Sw,So,Cp,aw,as,PR,Cpin,SGM,CMP,KWOG,QQ)

wmu=fusion(Cp(Won(:,1)),PR.mup);

TW=KWOG.w(Won(:,1))./CMP.Bw(Won(:,1),2)./wmu;
TO=KWOG.o(Won(:,1))./CMP.Bo(Won(:,1),2)/PR.mu(2);
TG=KWOG.g(Won(:,1))./CMP.Bg(Won(:,1),2)/PR.mu(3);
TP=KWOG.w(Won(:,1))./wmu.*Cp(Won(:,1));
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
Wg = Won(:,2).*TG.*(Uf==1);
W1 = SGM.Cg(Won(:,1)).*Wg + (CMP.Rs(Won(:,1),2).*SGM.Cor(Won(:,1)) + SGM.Co(Won(:,1))).*Wo + SGM.Cw(Won(:,1)).*W6;

W7=-Won(:,2).*SGM.Cw(Won(:,1)).*(TP.*(Uf==1)+Tip.*(Uf==-1));%

dP = dPw - dP(Won(:,1));
QQ = QQ + (SGM.Cg(Won(:,1)).*Wg + (CMP.Rs(Won(:,1),2).*SGM.Cor(Won(:,1)) + SGM.Co(Won(:,1))).*Wo + SGM.Cw(Won(:,1)).*W6).*dP;


