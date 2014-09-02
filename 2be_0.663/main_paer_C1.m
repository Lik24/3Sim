function SD=main_paer_C1

clear all
addpath('Sim_Lib','Tube_Lib','Gor_crack','Sparse_GPU','CrGeom','Termal_lib','GeoMeh_Lib','DATA');

PR=Gl_PRM;

%[KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z,XYc]=Sintetic_Real(PR.Ns,PR.Nl);
[KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z]=Sintetic(PR.Ns,PR.Nl);
[WData]=Well_DATA(WXY,Z,PR.Ta);
Sw(:)=0;
[nt,PXY,gXY,PR.dl,tXY]=kvad_crack_fun5(WXY,PR.Nl);
[DATA]=GridProp(KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z,gXY,PR.Nl,WXY);

GYData=GY_DATA(DATA.BndXY,DATA.BndZ);
[nt1,PXY1]=derevo(nt,DATA.XY,22);

[CrDATA]=CrackProp(DATA,PR.dl);

% [nt,PXY]=Tube_perc(PR,CrDATA,DATA.XY,1.1,WXY);

[gt,GS]=Tresh_Gor(1,DATA.XY,PR.Nl);

[C,A2C,dVc,pc,DATA.WonV]=Conek(DATA.XY,nt1,PR.Nl,CrDATA,DATA.Won,PR.dh,PR.Kc,WData.r0);

[nt2,PXY2]=derevo(nt,DATA.XY,24);

gt(:)={[]};
[G,A2G,dVg,pg,DATA.WonG]=Conek(DATA.XY,nt2,PR.Nl,CrDATA,DATA.Won,PR.dh,PR.Kc,WData.r0);%Gorizont(DATA.XY,GS,gt,WXY);

parfor i=1:5
   SD(i)=Par_Sim_MKT1(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,i,PXY1,PXY2);
end;
end


function CD=Par_Sim_MKT1(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,i,PXY1,PXY2)

[XY,KX,Z,Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,temp_c]=SimT_MKT(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,i);

SD(1,1)={Q};
SD(2,1)={Sw};
SD(3,1)={Pi};
SD(4,1)={XY};
%SD(5,1)={WXY};
SD(6,1)={MCp};
SD(7,1)={p};
SD(8,1)={PXY1};
SD(9,1)={temp_c};
SD(10,1)={Z};
SD(11,1)={NDT};
SD(12,1)={Pw};
SD(13,1)={PpW};

CD={SD};

we=Q(:,2,:);
i
sum(we(:))
end