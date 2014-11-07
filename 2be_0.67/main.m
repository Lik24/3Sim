clearvars

addpath('Sim_Lib','Tube_Lib','Gor_crack','Sparse_GPU','CrGeom','Termal_lib','GeoMeh_Lib',...
    'Well_lib','Crack_gen','Problems','Poly_lib','SS_lib','Diff_lib','Viz_lib','DATA_In','Adap_lib');
PR=Gl_PRM;

[KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z,XY_GY,XY_GY_new]=Sintetic_Real(PR.Ns,PR.Nl);
KX=10*KX;
KX(:)=1;
XY_GYs=XY_GY;
%[KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z]=Sintetic(PR.Ns,PR.Nl);
[WData]=Well_DATA(WXY,Z,PR.Ta);
%[WData,Ppwf,Pw_d,Pw_z]=Well_DATA_Adap(WXY,Z,PR.Ta);
Sw(:)=0;
[nt,PXY,gXY,PR.dl,tXY,XY_GY]=kvad_crack_fun(XY_GY,PR.Nl,WXY);
[DATA]=GridProp(KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY(1:33,:),H,Z,tXY,PR.Nl,WXY,XY_GY);
Sw0=DATA.gSw;
[GYData,DATA.gKX,DATA.gSw,B,A2B,dVb]=GY_DATA(DATA,XY_GYs,XY_GY_new,PR);
%[WData.Doly,DATA.gKX,GYData.GY_Kxy]=Load_adp_prm2(WData.Doly,DATA.gKX,GYData.GY_Kxy);
%[WData.Doly,DATA,GYData]=Load_adp_prm(DATA,GYData,tXY);
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

d2=[];
[D,A2D,dVd,pd,DATA.WonD,DATA.Ld]=DobPor(DATA.XY,d2,PR.Nl,DATA.Won,WData.r0,DATA.ka);

%[XY,K,Z,Pi,Sw,Cp,p,Q]=Sim_MKT(Prop,C,A2C,G,A2G,dVc,dVg,DATA,WData);
[Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,Uf,dt1,dV0]=SimT_MKT(PR,C,A2C,G,A2G,B,A2B,D,A2D,dVc,dVg,dVb,DATA,WData,GYData,1,CR_GRUP);

VZL(DATA,WXY,Pi,Sw(:,2),Ti,MCp,PR.Nl,p,Q,SwC,CR_GRUP,pc,nt,XY_GY,Uf(:,end));
%VZL_VORONOI(XY,Sw(:,end),p,WXY,WData.Uf(:,end))


  Qo(:,1)=sum(Q(:,3,:));
  Ql(:,1)=sum(Q(:,2,:));
  Qz(:,1)=sum(Q(:,1,:));

 sQo=cumsum(Qo,1);
 sQl=cumsum(Ql,1);


c=1-Qo./Ql;
dV0(p)=dV0;
V0=sum(dV0.*(1-Sw0));
sQo(end,:)/V0