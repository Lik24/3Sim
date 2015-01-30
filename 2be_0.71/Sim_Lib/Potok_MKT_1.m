function [TW,TP,TG]=Potok_MKT_1(T,vP,Sw,Cp,as,aw,mu,rc)

r=rc(:,1);
c=rc(:,2);

Kfw=Sat_cal(Sw,1,1,as,aw); %water

%Swc=Sw(c);
%Swl=Sw(r);

Kwc=Kfw(r);
Kwl=Kfw(c);

Cpc=Cp(r);
Cpl=Cp(c);

%Swe=Swc.*vP+Swl.*(vP==0);

Cpe=Cpc.*vP(:,1)+Cpl.*vP(:,2);
Kf=Kwc.*vP(:,1)+Kwl.*vP(:,2);

%Kf=Sat_cal(Swe,1,1,as,aw); %water
%Kfp=Sat_cal(Swe.*Cpe,1,1,as,aw); %polim

TW=T.*Kf.*((1-Cpe)./mu(1)+Cpe./mu(4));
TP=T.*Kf.*Cpe./mu(4);

TG=1;