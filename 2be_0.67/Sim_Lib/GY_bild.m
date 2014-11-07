function [bal,bcl,bgl,bdl,bbl,baw,bdw]=GY_bild(GY_Data,Pi,Sw,BXY,BZ,na,nc,ng,nd,nb,Txy_GY,Tz_GY,as,aw,mu)
    
GY_Pxy=GY_Data.GY_Pxy;
GY_Swxy=GY_Data.GY_Swxy;

GY_Pz=GY_Data.GY_Pz;
GY_Swz=GY_Data.GY_Swz;

GY_Txy=GY_Data.GY_Txy;



vP=GY_Pxy>=Pi;
SwT=GY_Swxy.*vP+Sw.*(vP==0);

Kfw_xy=Sat_cal(SwT,1,1,as,aw); %water
Kfo_xy=Sat_cal(SwT,2,1,as,aw); %oil

vP=GY_Pz>=Pi;
SwT=GY_Swz.*vP+Sw.*(vP==0);

Kfw_z=Sat_cal(SwT,1,1,as,aw); %water
Kfo_z=Sat_cal(SwT,2,1,as,aw); %oil

Tw_xy=Txy_GY.*Kfw_xy/mu(1);
To_xy=Txy_GY.*Kfo_xy/mu(2);

Twz=Tz_GY.*Kfw_z/mu(1);
Toz=Tz_GY.*Kfo_z/mu(2);

bal=zeros(2,na);
bal(2,:)=Twz+Toz;
bal(1,:)=Tw_xy'+To_xy';

baw=zeros(2,na);
baw(2,:)=Twz;
baw(1,:)=Tw_xy';

bcl=zeros(1,nc);
bgl=zeros(1,ng);
bdl=zeros(2,nd);
bdw=zeros(2,nd);


bbl=zeros(1,nb);
bbl(1,:)=GY_Txy.*aw(1)/mu(1);