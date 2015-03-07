function [W1,W6,Wo,Wg,W7]=Well_MKT(Won,Uf,Sw,Cp,aw,as,PR,Cpin,SGM)

wmu=fusion(Cp(Won(:,1)),PR.mup);
TW=Sat_cal(Sw(Won(:,1)),1,1,as,aw)./SGM.Bwog(Won(:,1),2)./wmu;
TO=Sat_cal(Sw(Won(:,1)),2,1,as,aw)./SGM.Bwog(Won(:,1),4)/PR.mu(2);
TG=Sat_cal(Sw(Won(:,1)),3,1,as,aw)./SGM.Bwog(Won(:,1),6)/PR.mu(3);
TP=Sat_cal(Sw(Won(:,1)),1,1,as,aw)./wmu.*Cp(Won(:,1));
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
W1 = Wg + (SGM.Rs(Won) - SGM.Cgso(Won)).*Wo - SGM.Cgsw(Won).*W6;%?????????????
W7=-Won(:,2).*SGM.Cgsw(Won).*(TP.*(Uf==1)+Tip.*(Uf==-1));%
