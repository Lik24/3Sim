function [TL,TW,TP,TG]=Potok_MKT_2(T,vP,Cp,mu,sd,na,kfw,kfo,kms,dP,L,Ro,Ke)

Kwl=kfw(sd(:,1)); %water
Kwc=kfw(sd(:,2)); %water

Kol=kfo(sd(:,1)); %lic
Koc=kfo(sd(:,2)); %lic

Cpl=Cp(sd(:,1));
Cpc=Cp(sd(:,2));

Cpe=Cpc.*vP(:,1)+Cpl.*vP(:,2);
Kfw=Kwc.*vP(:,1)+Kwl.*vP(:,2);
Kfo=Koc.*vP(:,1)+Kol.*vP(:,2);

TO=T.*Kfo./mu(2);
TW=T.*Kfw.*((1-Cpe)./mu(1)+Cpe./mu(4));
TP=T.*Kfw.*Cpe./mu(4);
TG=1;


 if kms~=0
    %     проверить првильность задания плотности
    dPL=abs(dP)./L(sd(:,1)+(sd(:,2)-1)*na);
    K=Ke(sd(:,1));
    [TW] = Forh(TW,kms, Ro(1), Kfw, K, mu(1), dPL);
    [TO] = Forh(TO,kms, Ro(2), Kfo, K, mu(2), dPL);
end;

TL=sparse(sd(:,1),sd(:,2),TW+TO,na,na);
TW=sparse(sd(:,1),sd(:,2),TW,na,na);

TL=TL+TL';
TW=TW+TW';

 