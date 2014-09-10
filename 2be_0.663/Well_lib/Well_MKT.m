function [W1,W6,W7]=Well_MKT(Wf,Won,Uf,Sw,Cp,aw,as,mu,mup,Cpin)

wmu=fusion(Cp(Won),mup);
TW=Sat_cal(Sw(Won),1,1,as,aw)./wmu;
TO=Sat_cal(Sw(Won),2,1,as,aw)/mu(2);
TP=Sat_cal(Sw(Won),1,1,as,aw)./wmu.*Cp(Won);
% 
% TW=kfw(Won)/mu(1);
% TO=kfo(Won),2,1,as,aw)/mu(2);
% TP=Sat_cal(Sw(Won),1,1,as,aw)/mu(4).*Cp(Won);

Tiw=(1-Cpin)./wmu;
Tip=Cpin./wmu;%Sat_cal(ones(size(Won,1),1).*Cpin,1,1,as,aw)/mu(4);

% size((TW+TO+TP))
% size((Uf==1))
% size((Tiw+Tip))
% size((Uf==-1))
%  size(Wf)
%  size(TW)
%  size(Uf)

W1=Wf.*((TW+TO).*(Uf==1)+(Tiw+Tip).*(Uf==-1));%
W6=Wf.*(TW.*(Uf==1)+(Tiw+Tip).*(Uf==-1));%
W7=Wf.*(TP.*(Uf==1)+Tip.*(Uf==-1));%
