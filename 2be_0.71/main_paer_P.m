function [SD,Q,Ca]=main_paer_P
addpath('Sim_Lib','Tube_Lib','Gor_crack','Sparse_GPU','CrGeom','Termal_lib','GeoMeh_Lib',...
    'Well_lib','Crack_gen','Problems','Poly_lib','SS_lib','Diff_lib','Viz_lib','DATA_In','Adap_lib');
PR=Gl_PRM;

[KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z,XYc]=Sintetic_Real(PR.Ns,PR.Nl);
KX=10*KX;
%[KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z]=Sintetic(PR.Ns,PR.Nl);
[WData]=Well_DATA(WXY,Z,PR.Ta);
%Sw(:)=0;
[nt,PXY,gXY,PR.dl,tXY,XYc]=kvad_crack_fun(XYc,PR.Nl,WXY);
[DATA]=GridProp(KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY(1:33,:),H,Z,tXY,PR.Nl,WXY,XYc);
Sw0=DATA.gSw;
GYData=GY_DATA(DATA.BndXY,DATA.BndZ,DATA.gP);
[WData.Doly,DATA,GYData]=Load_adp_prm(DATA,GYData,tXY);
%[nt1,PXY]=derevo(nt,DATA.XY,22);
%[nt,PXY]=elka(PR.Nl,DATA.XY,6,10,0,25);  % кол-во трещин, длинна, флаг к скважине

%load('elka_tst.mat','nt','PXY')
[CrDATA]=CrackProp(DATA,PR,nt);

% [nt,PXY]=Tube_perc(PR,CrDATA,DATA.XY,1.1,WXY);

[gt,GS]=Tresh_Gor(1,DATA.XY,PR.Nl);

nt0={[]};
nt(:)={nt0};
[C,A2C,dVc,pc,DATA.WonV,DATA.Lc,CR_GRUP]=Conek2(DATA.XY,nt,PR.Nl,CrDATA,DATA.Won,WData.r0,DATA.ka);

%[nt2,PXY2]=derevo(nt,DATA.XY,23);
gt(:)={[]};
%nt2(:)={[]};
[G,A2G,dVg,pg,DATA.WonG,DATA.Lg]=Conek(DATA.XY,gt,PR.Nl,CrDATA,DATA.Won,PR.dh,PR.Kc,WData.r0);%Gorizont(DATA.XY,GS,gt,WXY,WData.r0);

parfor i=1:3
   [SD(i),Q(:,i),Ca(:,i)]=Zap_exp(PR,C,A2C,G,A2G,dVc,dVg,DATA,GYData,CR_GRUP,Sw0,WXY,Z,i,XYc,tXY); 
end
%VZL(DATA,WXY,Pi,Sw,Ti,MCp,PR.Nl,p,Q,SwC,CR_GRUP,pc,nt,XYc,WData.Uf(:,end));
%VZL_VORONOI(XY,Sw(:,end),p,WXY,WData.Uf(:,end))
end
function [SD,Qaut,c]=Zap_exp(PR,C,A2C,G,A2G,dVc,dVg,DATA,GYData,CR_GRUP,Sw0,WXY,Z,i,XYc,tXY)
  [WData]=Well_DATA_mP(WXY,Z,PR.Ta,i);
  [WData.Doly,DATA,GYData]=Load_adp_prm(DATA,GYData,tXY);
  [Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,WData.Uf,dt1,dV0]=SimT_MKT(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,1,CR_GRUP);
  
  Qo(:,1)=sum(Q(:,3,:));
  Ql(:,1)=sum(Q(:,2,:));
  Qz(:,1)=sum(Q(:,1,:));

 sQo=cumsum(Qo,1);
 sQl=cumsum(Ql,1);

 Qaut=sQo;

c=1-Qo./Ql;
dV0(p)=dV0;
V0=sum(dV0.*(1-Sw0));
sQo(end,:)/V0

CD(1)={Q};
CD(2)={Pi};
CD(3)={Sw};
CD(4)={p};
CD(5)={PpW};
CD(6)={dV0};
CD(7)={DATA};
CD(8)={WXY};
CD(9)={WData};
CD(10)={Sw0};
CD(11)={PR};
CD(12)={XYc};
CD(13)={Pw};
SD={CD};
end