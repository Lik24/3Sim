clearvars

addpath('Sim_Lib','Tube_Lib','Gor_crack','Sparse_GPU','CrGeom','Termal_lib','GeoMeh_Lib',...
    'Well_lib','Crack_gen','Problems','Poly_lib','SS_lib','Diff_lib','Viz_lib','DATA_In','Adap_lib');
PR=Gl_PRM;

[KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z,XYc]=Sintetic_Real(PR.Ns,PR.Nl);
%[KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z]=Sintetic(PR.Ns,PR.Nl);
KX=KX*10;
[WData,Ppwf,Pw_d,Pw_z]=Well_DATA_Adap(WXY,Z,PR.Ta);
%Sw(:)=0;
[nt,PXY,gXY,PR.dl,tXY,XYc]=kvad_crack_fun(XYc,PR.Nl,WXY);
[DATA]=GridProp(KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z,tXY,PR.Nl,WXY,XYc);
GYData=GY_DATA(DATA.BndXY,DATA.BndZ,DATA.gP);
%[WData.Doly,DATA.gKX,GYData.GY_Kxy]=Load_adp_prm2(WData.Doly,DATA.gKX,GYData.GY_Kxy);

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
Sw0=Sw;
 Pw_f=Pw_d+Pw_z; 
%[Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,dt1,V0,Pw1,MJ,DATA.gKX,WData.Doly]=Adap_Kswaz_All_2(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,CR_GRUP,Pw_f);
AAS=load('ADP_ReZ3');
DATA.gKX=AAS.DATA.gKX;
WData.Doly=AAS.WData.Doly;
[Pi,Sw,Ti,MCp,p,Q,Pw4,PpW,SwC,NDT,dQ,dSS,dt1,V0,Pw3,MJ,A1,Kxy]=Adap_GY(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,CR_GRUP,Pw_f);


VZL(DATA,WXY,Pi,Sw,Ti,MCp,PR.Nl,p,Q,SwC,CR_GRUP,pc,nt,XYc,WData.Uf(:,end));
%VZL_VORONOI(XY,Sw(:,end),p,WXY,WData.Uf(:,end))


  Qo(:,1)=sum(Q(:,3,:));
  Ql(:,1)=sum(Q(:,2,:));
  Qz(:,1)=sum(Q(:,1,:));

 sQo=cumsum(Qo,1);
 sQl=cumsum(Ql,1);


c=1-Qo./Ql;
sQo(end,:)/sum(V0)