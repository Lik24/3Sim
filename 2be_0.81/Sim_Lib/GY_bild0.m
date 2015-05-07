function [WBND,QQBND,QQwoBND,KWOG_GY]=GY_bild(GY_Data,Pi,Sw,So,Cp,RC,T_GY_A,T_GY_D,PR,SGM,CMP,BXYZ,V)

as=PR.as;
aw=PR.aw;
ts=PR.ts;
tw=PR.tw;
mu=PR.mu;

vad=RC.ADr;

MSw=Sw(V.va);
DSw=Sw(V.vd);

MSo=So(V.va);
DSo=So(V.vd);

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
SoT = (1-GY_Swxy).*vP+MSo.*(vP==0);
MCpT = 0.*vP+MCp.*(vP==0);

KWOG_GY.Kfw_xyM=Sat_cal(SwT,SoT,1,1,as,aw); %water
KWOG_GY.Kfo_xyM=Sat_cal(SwT,SoT,2,1,as,aw); %oil
KWOG_GY.Kfg_xyM=Sat_cal(SwT,SoT,3,1,as,aw); %gas

vP=GY_Pz>=Pi(V.va);
SwT=GY_Swz.*vP+MSw.*(vP==0);

KWOG_GY.Kfw_zM=Sat_cal(SwT,SoT,1,1,as,aw); %water
KWOG_GY.Kfo_zM=Sat_cal(SwT,SoT,2,1,as,aw); %oil
KWOG_GY.Kfg_zM=Sat_cal(SwT,SoT,3,1,as,aw); %gas

%% Для двойной среды
vP=GY_Pxy(vad)>=Pi(V.vd);
SwT=GY_Swxy(vad).*vP+DSw.*(vP==0);
SoT = DSo.*(vP==0);
DCpT=0.*vP+DCp.*(vP==0);

KWOG_GY.Kfw_xyD=Sat_cal(SwT,SoT,1,1,ts,tw); %water
KWOG_GY.Kfo_xyD=Sat_cal(SwT,SoT,2,1,ts,tw); %oil
KWOG_GY.Kfg_xyD=Sat_cal(SwT,SoT,3,1,ts,tw); %gas

vP=GY_Pz(vad)>=Pi(V.vd);
SwT=GY_Swz(vad).*vP+DSw.*(vP==0);

KWOG_GY.Kfw_zD=Sat_cal(SwT,SoT,1,1,ts,tw); %water
KWOG_GY.Kfo_zD=Sat_cal(SwT,SoT,2,1,ts,tw); %oil
KWOG_GY.Kfg_zD=Sat_cal(SwT,SoT,3,1,ts,tw); %gas
%%
WBND.b1gm=Abra_GY(T_GY_A,KWOG_GY.Kfw_xyM,KWOG_GY.Kfo_xyM,KWOG_GY.Kfg_xyM,KWOG_GY.Kfw_zM,KWOG_GY.Kfo_zM,KWOG_GY.Kfg_zM,mu,MCpT,SGM,CMP,V.va);
WBND.b1gd=Abra_GY(T_GY_D,KWOG_GY.Kfw_xyD,KWOG_GY.Kfo_xyD,KWOG_GY.Kfg_xyD,KWOG_GY.Kfw_zD,KWOG_GY.Kfo_zD,KWOG_GY.Kfg_zD,mu,DCpT,SGM,CMP,V.vd);

WBND.b1gb = zeros(size(V.vb,1),1);
WBND.b1gb(V.vb,1) = GY_Txy(V.vb).*aw(1)/mu(1);

QQBND.Qm = WBND.b1gm(:,1).*(GY_Pxy - Pi(V.va)) + WBND.b1gm(:,2).*(GY_Pz - Pi(V.va));
QQBND.Qd = WBND.b1gd(:,1).*(GY_Pxy(vad) - Pi(V.vd)) + WBND.b1gd(:,2).*(GY_Pz(vad) - Pi(V.vd));
QQBND.Qb = WBND.b1gb.*(GY_Data.P0(V.vb) - Pi(V.vb))+ WBND.b1gb(V.vb,1).*(GY_Data.P0(V.vb) - Pi(V.vb));

QQwoBND.Qmw = WBND.b1gm(:,3).*(GY_Pxy - Pi(V.va)) + WBND.b1gm(:,4).*(GY_Pz - Pi(V.va));
QQwoBND.Qmo = WBND.b1gm(:,5).*(GY_Pxy - Pi(V.va)) + WBND.b1gm(:,6).*(GY_Pz - Pi(V.va));
QQwoBND.Qdw = WBND.b1gd(:,3).*(GY_Pxy(vad) - Pi(V.vd)) + WBND.b1gd(:,4).*(GY_Pz(vad) - Pi(V.vd));
QQwoBND.Qdo = WBND.b1gd(:,5).*(GY_Pxy(vad) - Pi(V.vd)) + WBND.b1gd(:,6).*(GY_Pz(vad) - Pi(V.vd)); 
end

function bal=Abra_GY(T_GY,Kfw_xy,Kfo_xy,Kfg_xy,Kfw_z,Kfo_z,Kfg_z,mu,Cp,SGM,CMP,v)

Tw_xy=T_GY(:,1).*Kfw_xy/mu(1)./CMP.Bw(v,2);
To_xy=T_GY(:,1).*Kfo_xy/mu(2)./CMP.Bo(v,2);
Tg_xy=T_GY(:,1).*Kfg_xy/mu(3)./CMP.Bg(v,2);

Twz=T_GY(:,2).*Kfw_z/mu(1)./CMP.Bw(v,2);
Toz=T_GY(:,2).*Kfo_z/mu(2)./CMP.Bo(v,2);
Tgz=T_GY(:,2).*Kfo_z/mu(3)./CMP.Bg(v,2);

bal(:,1) = Tg_xy + (CMP.Rs(v,2).*SGM.Cor(v) + SGM.Co(v)).*To_xy + SGM.Cw(v).*Tw_xy;
bal(:,2) = Tgz + (CMP.Rs(v,2).*SGM.Cor(v) + SGM.Co(v)).*Toz + SGM.Cw(v).*Twz;
bal(:,3)=Tw_xy;
bal(:,4)=Twz;
bal(:,5)=To_xy;
bal(:,6)=Toz;
bal(:,7)=Tw_xy.*Cp;
bal(:,8)=Twz.*Cp;
end