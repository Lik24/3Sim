function [W1,W6,W7]=Well_MKT_2(Wf,Won,Uf,Cp,mu,Cpin,kfw,kfo)

% TW=Sat_cal(Sw(Won),1,1,as,aw)/mu(1);
% TO=Sat_cal(Sw(Won),2,1,as,aw)/mu(2);
% TP=Sat_cal(Sw(Won),1,1,as,aw)/mu(4).*Cp(Won);

TW=kfw(Won)/mu(1);
TO=kfo(Won)/mu(2);
TP=kfw(Won)/mu(4).*Cp(Won);

Tiw=(1-Cpin)/mu(1);%Sat_cal(ones(size(Won,1),1).*(1-Cpin),1,1,as,aw)/mu(1);
Tip=Cpin/mu(4);%Sat_cal(ones(size(Won,1),1).*Cpin,1,1,as,aw)/mu(4);

% size((TW+TO+TP))
% size(Wf)
% size(TW+TO)
% size(Uf)

W1=Wf.*((TW+TO).*(Uf==1)+(Tiw+Tip).*(Uf==-1));%
W6=Wf.*(TW.*(Uf==1)+(Tiw+Tip).*(Uf==-1));%
W7=Wf.*(TP.*(Uf==1)+Tip.*(Uf==-1));%
