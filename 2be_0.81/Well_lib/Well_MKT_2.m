function [W1,W6,Wo,W7]=Well_MKTF2(Wf,Won,Uf,Cp,mu,Cpin,kfw,kfo)

TW=kfw(Won)./SGM.Bwog(Won,2)/mu(1);
TO=kfo(Won)./SGM.Bwog(Won,4)/mu(2);
TP=kfw(Won)/mu(4).*Cp(Won);

Tiw=(1-Cpin)/mu(1);%Sat_cal(ones(size(Won,1),1).*(1-Cpin),1,1,as,aw)/mu(1);
Tip=Cpin/mu(4);%Sat_cal(ones(size(Won,1),1).*Cpin,1,1,as,aw)/mu(4);

% size((TW+TO+TP))
% size(Wf)
% size(TW+TO)
% size(Uf)


Wo = Wf.*TO.*(Uf==1);
W6 = Wf.*(TW.*(Uf==1)+(Tiw+Tip).*(Uf==-1));
W1 = SGM.Cg(Won).*Wf.*TG.*(Uf==1) + (SGM.Rs(Won).*SGM.Cor(Won) + SGM.Co(Won)).*Wo + SGM.Cw(Won).*W6
W7=Wf.*(TP.*(Uf==1)+Tip.*(Uf==-1));%
