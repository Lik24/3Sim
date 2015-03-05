function [TL,TW,TO,TG,TWCgs,TOCgs,TP]=Potok_MKT_2(T,vP,Cp,mu,sd,na,kfw,kfo,kfg,kms,dP,L,Ro,Ke,SGM,v1)

Kwl=kfw(sd(:,1)); %water
Kwc=kfw(sd(:,2)); %water

Kol=kfo(sd(:,1)); %lic
Koc=kfo(sd(:,2)); %lic

Kgl=kfg(sd(:,1)); %lic
Kgc=kfg(sd(:,2)); %lic

Cpl=Cp(sd(:,1));
Cpc=Cp(sd(:,2));

Cpe = Cpc.*vP(:,1) + Cpl.*vP(:,2);
Kfw = Kwc.*vP(:,1) + Kwl.*vP(:,2);
Kfo = Koc.*vP(:,3) + Kol.*vP(:,4);
Kfg = Kgc.*vP(:,5) + Kgl.*vP(:,6); 

Bw = 0.5*(SGM.Bwog(sd(:,1),2) + SGM.Bwog(sd(:,2),2));
Bo = 0.5*(SGM.Bwog(sd(:,1),4) + SGM.Bwog(sd(:,2),4));
Bg = 0.5*(SGM.Bwog(sd(:,1),6) + SGM.Bwog(sd(:,2),6));
Rs = 0.5*(SGM.Rs(sd(:,1)) + SGM.Rs(sd(:,2)));

TO = T.*Kfo./Bo./mu(2);
TW = T.*Kfw./Bw.*((1-Cpe)./mu(1)+Cpe./mu(4));
TP = T.*Kfw.*Cpe./mu(4);
TG = T.*Kfg./Bg./mu(3);

 if kms~=0
    %     проверить првильность задания плотности
    dPL=abs(dP)./L(sd(:,1)+(sd(:,2)-1)*na);
    K=Ke(sd(:,1));
    [TW] = Forh(TW,kms, Ro(1), Kfw, K, mu(1), dPL);
    [TO] = Forh(TO,kms, Ro(2), Kfo, K, mu(2), dPL);
end;

TW=sparse(sd(:,1),sd(:,2),TW,na,na);  TW=TW+TW';
TO=sparse(sd(:,1),sd(:,2),TO,na,na);  TO=TO+TO';
TG=sparse(sd(:,1),sd(:,2),TG,na,na);  TG=TG+TG';

TW=TW(:,v1==1);
TW=TW(v1==1,:);

TO=TO(:,v1==1);
TO=TO(v1==1,:);

TG=TG(:,v1==1);
TG=TG(v1==1,:);
 
n = size(TO,1);
[r,c]=find(TO+TW);

TWa = TW.*sparse(r,c,-SGM.Cgsw(c));
TOa = TO.*sparse(r,c,Rs(c)-SGM.Cgso(c));

TW = TW - sparse(1:n,1:n,sum(TW,2),n,n);
TO = TO - sparse(1:n,1:n,sum(TO,2),n,n);

TWCgs = TWa - sparse(1:n,1:n,sum(TWa,2),n,n);
TOCgs = TOa - sparse(1:n,1:n,sum(TOa,2),n,n);
TG = TG - sparse(1:n,1:n,sum(TG,2),n,n);
TL = TG + TWCgs + TOCgs;
