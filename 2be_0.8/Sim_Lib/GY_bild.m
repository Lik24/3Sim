function [BOUNDFL]=GY_bild(GY_Data,Pi,Sw,Cp,RC,T_GY_A,T_GY_D,PR)

as=PR.as;
aw=PR.aw;
ts=PR.ts;
tw=PR.tw;
mu=PR.mu;

na=RC.na;   nc=RC.nc;   ng=RC.ng;   nd=RC.nd;  nb=RC.nb;  
va=1:na;
vd=na+1:na+nd;
vad=RC.ADr;

MSw=Sw(1:na);
DSw=Sw(na+1:end);

MCp=Cp(1:na);
DCp=Cp(na+1:end);

GY_Pxy=GY_Data.GY_Pxy;
GY_Swxy=GY_Data.GY_Swxy;

GY_Pz=GY_Data.GY_Pz;
GY_Swz=GY_Data.GY_Swz;

GY_Txy=GY_Data.GY_Txy;
%% Для матрицы
vP=GY_Pxy(va)>=Pi(va);
SwT=GY_Swxy.*vP+MSw.*(vP==0);
MCpT=0.*vP+MCp.*(vP==0);

Kfw_xyM=Sat_cal(SwT,1,1,as,aw); %water
Kfo_xyM=Sat_cal(SwT,2,1,as,aw); %oil

vP=GY_Pz>=Pi(1:na);
SwT=GY_Swz.*vP+MSw.*(vP==0);

Kfw_zM=Sat_cal(SwT,1,1,as,aw); %water
Kfo_zM=Sat_cal(SwT,2,1,as,aw); %oil

%% Для двойной среды
vP=GY_Pxy(vad)>=Pi(vd);
SwT=GY_Swxy(vad).*vP+DSw.*(vP==0);
DCpT=0.*vP+DCp.*(vP==0);

Kfw_xyD=Sat_cal(SwT,1,1,ts,tw); %water
Kfo_xyD=Sat_cal(SwT,2,1,ts,tw); %oil

vP=GY_Pz(vad)>=Pi(vd);
SwT=GY_Swz(vad).*vP+DSw.*(vP==0);

Kfw_zD=Sat_cal(SwT,1,1,ts,tw); %water
Kfo_zD=Sat_cal(SwT,2,1,ts,tw); %oil
%%
BOUNDFL.b1gm=Abra_GY(T_GY_A,Kfw_xyM,Kfo_xyM,Kfw_zM,Kfo_zM,mu,MCpT);
BOUNDFL.b1gd=Abra_GY(T_GY_D,Kfw_xyD,Kfo_xyD,Kfw_zD,Kfo_zD,mu,DCpT);

BOUNDFL.b1gc=zeros(nc,1);
BOUNDFL.b1gg=zeros(ng,1);

BOUNDFL.b1gb=zeros(nb,1);
BOUNDFL.b1gb(:,1)=GY_Txy.*aw(1)/mu(1);
end

function bal=Abra_GY(T_GY,Kfw_xy,Kfo_xy,Kfw_z,Kfo_z,mu,Cp)

Tw_xy=T_GY(:,1).*Kfw_xy/mu(1);
To_xy=T_GY(:,1).*Kfo_xy/mu(2);

Twz=T_GY(:,2).*Kfw_z/mu(1);
Toz=T_GY(:,2).*Kfo_z/mu(2);

bal(:,1)=Tw_xy+To_xy;
bal(:,2)=Twz+Toz;
bal(:,3)=Tw_xy;
bal(:,4)=Twz;
bal(:,5)=Tw_xy.*Cp;
bal(:,6)=Twz.*Cp;
end