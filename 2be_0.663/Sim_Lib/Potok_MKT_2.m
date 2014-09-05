function [TL,TW,TP,TG]=Potok_MKT_2(T,vP,Cp,mu,sd,na,kfw,kfo)

Kwl=kfw(sd(:,1)); %water
Kwc=kfw(sd(:,2)); %water

Kol=kfo(sd(:,1)); %lic
Koc=kfo(sd(:,2)); %lic

Cpl=Cp(sd(:,1));
Cpc=Cp(sd(:,2));

Cpe=Cpc.*vP(:,1)+Cpl.*vP(:,2);
Kf=Kwc.*vP(:,1)+Kwl.*vP(:,2);
Kfo=Koc.*vP(:,1)+Kol.*vP(:,2);

TO=T.*Kfo./mu(2);
TW=T.*Kf.*((1-Cpe)./mu(1)+Cpe./mu(4));
TP=T.*Kf.*Cpe./mu(4);

TG=1;

TL=sparse(sd(:,1),sd(:,2),TW+TO,na,na);
TW=sparse(sd(:,1),sd(:,2),TW,na,na);

TL=TL+TL';
TW=TW+TW';

 