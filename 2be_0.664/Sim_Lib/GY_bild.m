function [bal,bcl,bgl,baw]=GY_bild(GY_Data,Pi,Sw,BZ,na,nc,ng,T_GY,as,aw,mu)
    
GY_Pz=GY_Data.GY_Pz;
GY_Swz=GY_Data.GY_Swz;

vP=GY_Pz(BZ~=0)>=Pi(BZ~=0);
SwT=GY_Swz(BZ~=0).*vP+Sw(BZ~=0).*(vP==0);

Kfw=Sat_cal(SwT,1,1,as,aw); %water
Kfo=Sat_cal(SwT,2,1,as,aw); %oil


Tw=T_GY(BZ~=0).*Kfw/mu(1);
To=T_GY(BZ~=0).*Kfo/mu(2);

bal=zeros(1,na);
bal(BZ~=0)=Tw+To;

baw=zeros(1,na);
baw(BZ~=0)=Tw;

GY_Kxy=GY_Data.GY_Kxy;
GY_Pxy=GY_Data.GY_Pxy;
GY_Swxy=GY_Data.GY_Swxy;


bcl=zeros(1,nc);
bgl=zeros(1,ng);