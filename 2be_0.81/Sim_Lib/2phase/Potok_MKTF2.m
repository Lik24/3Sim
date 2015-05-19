function [TW,TO,TP]=Potok_MKTF2(T,P,Cp,mu,sd,na,kfw,kfo,kms,L,Ro,Ke,CMP,v)

PW = P(v,1);
PO = P(v,2);
dPW = PW(sd(:,2)) - PW(sd(:,1));

vPW = dPW>=0; 
vPO = PO(sd(:,2)) - PO(sd(:,1))>=0;

Kwl=kfw(sd(:,1)); %water
Kwc=kfw(sd(:,2)); %water

Kol=kfo(sd(:,1)); %lic
Koc=kfo(sd(:,2)); %lic

Cpl=Cp(sd(:,1));
Cpc=Cp(sd(:,2));

Cpe = Cpc.*vPW + Cpl.*(vPW==0);
Kfw = Kwc.*vPW + Kwl.*(vPW==0);
Kfo = Koc.*vPO + Kol.*(vPO==0); 

Bw = 0.5*(CMP.Bw(v(sd(:,1)),2) + CMP.Bw(v(sd(:,2)),2));
Bo = 0.5*(CMP.Bo(v(sd(:,1)),2) + CMP.Bo(v(sd(:,2)),2));

TO = T.*Kfo./Bo./mu(2);
TW = T.*Kfw./Bw.*((1-Cpe)./mu(1)+Cpe./mu(4));
TP = T.*Kfw.*Cpe./mu(4);

 if kms~=0
    %     проверить првильность задания плотности
    dPL=abs(dPW)./L(sd(:,1)+(sd(:,2)-1)*na);
    K=Ke(sd(:,1));
    [TW] = Forh(TW,kms, Ro(1), Kfw, K, mu(1), dPL);
    [TO] = Forh(TO,kms, Ro(2), Kfo, K, mu(2), dPL);
end;

TW=sparse(sd(:,1),sd(:,2),TW,na,na);  TW=TW+TW';
TO=sparse(sd(:,1),sd(:,2),TO,na,na);  TO=TO+TO';

TW = TW - sparse(1:na,1:na,sum(TW,2),na,na);
TO = TO - sparse(1:na,1:na,sum(TO,2),na,na);

