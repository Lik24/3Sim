function [W1,W6,W7]=Well_MKT(Wf,Won,Uf,Sw,Cp,aw,as,mu,mup,Cpin,A,Bwo)

wmu=fusion(Cp(Won),mup);
TW=Sat_cal(Sw(Won),1,1,as,aw)./wmu./Bwo(Won,1);
TO=Sat_cal(Sw(Won),2,1,as,aw)/mu(2)./Bwo(Won,2);
TP=Sat_cal(Sw(Won),1,1,as,aw)./wmu.*Cp(Won);
% 
% TW=kfw(Won)/mu(1);
% TO=kfo(Won),2,1,as,aw)/mu(2);
% TP=Sat_cal(Sw(Won),1,1,as,aw)/mu(4).*Cp(Won);

Tiw=(1-Cpin)./wmu./Bwo(Won,1);
Tip=Cpin./wmu;

% size((TW+TO+TP))
% size((Uf==1))
% size((Tiw+Tip))
% size((Uf==-1))
% size(Wf)
% size(TW)
% size(Uf)
%A(:)=1;
W1=Wf.*((A(Won).*TW+TO).*(Uf==1)+A(Won).*(Tiw+Tip).*(Uf==-1));%
W6=Wf.*(TW.*(Uf==1)+(Tiw+Tip).*(Uf==-1));%./Bwo(Won,2)
W7=Wf.*A(Won).*(TP.*(Uf==1)+Tip.*(Uf==-1));%
