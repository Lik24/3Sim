clearvars

addpath('Sim_Lib','Tube_Lib','Gor_crack','Sparse_GPU','CrGeom','Termal_lib','GeoMeh_Lib',...
    'DATA','Well_lib','Crack_gen','Problems','Poly_lib','SS_lib','Diff_lib','Viz_lib');
PR=Gl_PRM;

%[KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z,XYc]=Sintetic_Real(PR.Ns,PR.Nl);
[KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z]=Sintetic(PR.Ns,PR.Nl);
[WData]=Well_DATA(WXY,Z,PR.Ta);
%Sw(:)=0;
[nt,PXY,gXY,PR.dl,tXY]=kvad_crack_fun(WXY,PR.Nl);
[DATA]=GridProp(KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z,gXY,PR.Nl,WXY);

GYData=GY_DATA(DATA.BndXY,DATA.BndZ);
%[nt1,PXY]=derevo(nt,DATA.XY,22);
%[nt,PXY]=elka(PR.Nl,DATA.XY,10,20,0,25);  % кол-во трещин, длинна, флаг к скважине

%load('elka_tst.mat','nt','PXY')
[CrDATA]=CrackProp(DATA,PR,nt);

% [nt,PXY]=Tube_perc(PR,CrDATA,DATA.XY,1.1,WXY);

[gt,GS]=Tresh_Gor(1,DATA.XY,PR.Nl);

%nt0={[]};
%nt(:)={nt0};
[C,A2C,dVc,pc,DATA.WonV,DATA.Lc,CR_GRUP]=Conek2(DATA.XY,nt,PR.Nl,CrDATA,DATA.Won,WData.r0);

%[nt2,PXY2]=derevo(nt,DATA.XY,23);
gt(:)={[]};
%nt2(:)={[]};
[G,A2G,dVg,pg,DATA.WonG,DATA.Lg]=Conek(DATA.XY,gt,PR.Nl,CrDATA,DATA.Won,PR.dh,PR.Kc,WData.r0);%Gorizont(DATA.XY,GS,gt,WXY,WData.r0);


%[XY,K,Z,Pi,Sw,Cp,p,Q]=Sim_MKT(Prop,C,A2C,G,A2G,dVc,dVg,DATA,WData);
[XY,KX,Z,Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,dt1,V0]=SimT_MKT(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,1,CR_GRUP);

VZL(XY,KX,WXY,Z,Pi,Sw,Ti,MCp,PR.Nl,p,Q,SwC,CR_GRUP,pc,nt);
%VZL_VORONOI(XY,Sw(:,end),p,WXY,WData.Uf(:,end))


  Qo(:,1)=sum(Q(:,3,:));
  Ql(:,1)=sum(Q(:,2,:));
  Qz(:,1)=sum(Q(:,1,:));

 sQo=cumsum(Qo,1);
 sQl=cumsum(Ql,1);


c=1-Qo./Ql;
sQo(end,:)/V0