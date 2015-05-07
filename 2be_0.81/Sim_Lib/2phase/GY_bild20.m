function [WBND,QQBND,QQwoBND,KWOG_GY]=GY_bild20(GY_Data,Pi,Sw,Cp,RC,T_GY_A,T_GY_D,PR,CMP,BXYZ,V)

as=PR.as;
aw=PR.aw;
ts=PR.ts;
tw=PR.tw;
mu=PR.mu;

vad=RC.ADr;

MSw=Sw(V.va);
DSw=Sw(V.vd);

MCp=Cp(V.va);
DCp=Cp(V.vd);

GY_Pxy=GY_Data.GY_Pxy;
GY_Swxy=GY_Data.GY_Swxy;

GY_Pz=GY_Data.GY_Pz;
GY_Swz=GY_Data.GY_Swz;

GY_Txy=GY_Data.GY_Txy;
%% Для матрицы
vP = GY_Pxy(V.va)>=Pi(V.va);
SwT = GY_Swxy.*vP+MSw.*(vP==0);
MCpT = 0.*vP+MCp.*(vP==0);

KWOG_GY.Kfw_xyM=Sat_cal(SwT,1 - SwT,1,1,as,aw); %water
KWOG_GY.Kfo_xyM=Sat_cal(SwT,1 - SwT,2,1,as,aw); %oil

vP=GY_Pz>=Pi(V.va);
SwT=GY_Swz.*vP+MSw.*(vP==0);

KWOG_GY.Kfw_zM=Sat_cal(SwT,1 - SwT,1,1,as,aw); %water
KWOG_GY.Kfo_zM=Sat_cal(SwT,1 - SwT,2,1,as,aw); %oil

%% Для двойной среды
vP=GY_Pxy(vad)>=Pi(V.vd);
SwT=GY_Swxy(vad).*vP+DSw.*(vP==0);
DCpT=0.*vP+DCp.*(vP==0);

KWOG_GY.Kfw_xyD=Sat_cal(SwT,1-SwT,1,1,ts,tw); %water
KWOG_GY.Kfo_xyD=Sat_cal(SwT,1-SwT,2,1,ts,tw); %oil

vP=GY_Pz(vad)>=Pi(V.vd);
SwT=GY_Swz(vad).*vP+DSw.*(vP==0);

KWOG_GY.Kfw_zD=Sat_cal(SwT,1-SwT,1,1,ts,tw); %water
KWOG_GY.Kfo_zD=Sat_cal(SwT,1-SwT,2,1,ts,tw); %oil

%%
WBND.b1gm=Abra_GY(T_GY_A,KWOG_GY.Kfw_xyM,KWOG_GY.Kfo_xyM,KWOG_GY.Kfw_zM,KWOG_GY.Kfo_zM,mu,MCpT,CMP,V.va);
WBND.b1gd=Abra_GY(T_GY_D,KWOG_GY.Kfw_xyD,KWOG_GY.Kfo_xyD,KWOG_GY.Kfw_zD,KWOG_GY.Kfo_zD,mu,DCpT,CMP,V.vd);

WBND.b1gb = zeros(size(V.vb,1),1);
WBND.b1gb(V.vb,1) = GY_Txy(V.vb).*aw(1)/mu(1);

QQBND.Qm = WBND.b1gm(:,1).*(GY_Pxy - Pi(V.va)) + WBND.b1gm(:,2).*(GY_Pz - Pi(V.va));
QQBND.Qd = WBND.b1gd(:,1).*(GY_Pxy(vad) - Pi(V.vd)) + WBND.b1gd(:,2).*(GY_Pz(vad) - Pi(V.vd));
QQBND.Qb = WBND.b1gb.*(GY_Data.P0(V.vb) - Pi(V.vb))+ WBND.b1gb(V.vb,1).*(GY_Data.P0(V.vb) - Pi(V.vb));

QQwoBND.Qmw = WBND.b1gm(:,3).*(GY_Pxy - Pi(V.va)) + WBND.b1gm(:,4).*(GY_Pz - Pi(V.va));
QQwoBND.Qdw = WBND.b1gd(:,3).*(GY_Pxy(vad) - Pi(V.vd)) + WBND.b1gd(:,4).*(GY_Pz(vad) - Pi(V.vd)); 
end

function bal=Abra_GY(T_GY,Kfw_xy,Kfo_xy,Kfw_z,Kfo_z,mu,Cp,CMP,v)

Tw_xy=T_GY(:,1).*Kfw_xy/mu(1)./CMP.Bw(v,2);
To_xy=T_GY(:,1).*Kfo_xy/mu(2)./CMP.Bo(v,2);

Twz=T_GY(:,2).*Kfw_z/mu(1)./CMP.Bw(v,2);
Toz=T_GY(:,2).*Kfo_z/mu(2)./CMP.Bo(v,2);

bal(:,1) = CMP.Cw(v).*Tw_xy + To_xy;
bal(:,2) = CMP.Cw(v).*Twz + Toz;
bal(:,3)=Tw_xy;
bal(:,4)=Twz;
bal(:,5)=To_xy;
bal(:,6)=Toz;
bal(:,7)=Tw_xy.*Cp;
bal(:,8)=Twz.*Cp;
end