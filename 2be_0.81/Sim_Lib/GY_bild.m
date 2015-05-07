function [WBND,QQBND]=GY_bild(GY_Data,Pi,dPt,KWOG_GY,Cp,RC,T_GY_A,T_GY_D,PR,SGM,CMP,V,QQBND)

mu=PR.mu;
MCp=Cp(V.va);
DCp=Cp(V.vd);

vad=RC.ADr;

GY_Pxy=GY_Data.GY_Pxy;
GY_Swxy=GY_Data.GY_Swxy;

GY_Pz=GY_Data.GY_Pz;
GY_Swz=GY_Data.GY_Swz;

GY_Txy=GY_Data.GY_Txy;

vP=GY_Pxy(V.va)>=Pi(V.va);
MCpT=0.*vP+MCp.*(vP==0);

vP=GY_Pxy(vad)>=Pi(V.vd);
DCpT=0.*vP+DCp.*(vP==0);

WBND.b1gm=Abra_GY(T_GY_A,KWOG_GY.Kfw_xyM,KWOG_GY.Kfo_xyM,KWOG_GY.Kfg_xyM,KWOG_GY.Kfw_zM,KWOG_GY.Kfo_zM,KWOG_GY.Kfg_zM,mu,MCpT,SGM,CMP,V.va);
WBND.b1gd=Abra_GY(T_GY_D,KWOG_GY.Kfw_xyD,KWOG_GY.Kfo_xyD,KWOG_GY.Kfg_xyD,KWOG_GY.Kfw_zD,KWOG_GY.Kfo_zD,KWOG_GY.Kfg_zD,mu,DCpT,SGM,CMP,V.vd);

WBND.b1gb = zeros(size(V.vb,1),1);
WBND.b1gb(V.vb,1) = GY_Txy(V.vb).*PR.aw(1)/mu(1);

QQBND.Qm = QQBND.Qm - WBND.b1gm(:,1).*dPt(V.va) - WBND.b1gm(:,2).*dPt(V.va);
QQBND.Qd = QQBND.Qd - WBND.b1gd(:,1).*dPt(V.vd) - WBND.b1gd(:,2).*dPt(V.vd);
QQBND.Qb = QQBND.Qb - WBND.b1gb.*dPt(V.vb) - WBND.b1gb(V.vb,1).*dPt(V.vb);

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