clear all
addpath('Sim_Lib','Tube_Lib','Gor_crack','Sparse_GPU','CrGeom','Termal_lib','GeoMeh_Lib',...
    'DATA','Well_lib','Crack_gen','Problems');
PR=Gl_PRM;

%[KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z,XYc]=Sintetic_Real(PR.Ns,PR.Nl);
[KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z]=Sintetic(PR.Ns,PR.Nl);
[WData]=Well_DATA(WXY,Z,PR.Ta);
Sw(:)=0;
[nt,PXY,gXY,PR.dl,tXY]=kvad_crack_fun5(WXY,PR.Nl);
[DATA]=GridProp(KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z,gXY,PR.Nl,WXY);

GYData=GY_DATA(DATA.BndXY,DATA.BndZ);
%[nt1,PXY]=derevo(nt,DATA.XY,22);
[nt1,PXY]=elka(nt,DATA.XY,10,10,0,25);  % ���-�� ������, ������, ���� � ��������

[CrDATA]=CrackProp(DATA,PR.dl);

% [nt,PXY]=Tube_perc(PR,CrDATA,DATA.XY,1.1,WXY);

[gt,GS]=Tresh_Gor(1,DATA.XY,PR.Nl);
%nt(:)={[]};
[C,A2C,dVc,pc,DATA.WonV]=Conek(DATA.XY,nt1,PR.Nl,CrDATA,DATA.Won,PR.dh,PR.Kc,WData.r0);

%[nt2,PXY2]=derevo(nt,DATA.XY,23);
gt(:)={[]};
%nt2(:)={[]};
[G,A2G,dVg,pg,DATA.WonG]=Gorizont(DATA.XY,GS,gt,WXY);%Conek(DATA.XY,nt2,PR.Nl,CrDATA,DATA.Won,PR.dh,PR.Kc,WData.r0);


%[XY,K,Z,Pi,Sw,Cp,p,Q]=Sim_MKT(Prop,C,A2C,G,A2G,dVc,dVg,DATA,WData);
[XY,KX,Z,Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS]=SimT_MKT(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,1);

VZL(XY,KX,WXY,Z,Pi,Sw,Ti,MCp,PR.Nl,p,Q,PXY,gt);

Qo=sum(Q(:,3,:));
sum(Qo)/(250*250*10*0.2)