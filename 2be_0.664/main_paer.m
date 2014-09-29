function SD=main_paer
clear all
addpath('Sim_Lib','Tube_Lib','Gor_crack','Sparse_GPU','CrGeom','Termal_lib','GeoMeh_Lib','DATA');

PR=Gl_PRM;

[KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z,XYc]=Sintetic_Real(PR.Ns,PR.Nl);

[nt,PXY,gXY,PR.dl]=kvad_crack_fun5(WXY,PR.Nl);

%nt(:)={[]};
%Sw(:)=0;
[DATA]=GridProp(KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,XYc,H,Z,gXY,PR.Nl,WXY);

 GYData=GY_DATA(DATA.BndXY,DATA.BndZ);
%[nt1,PXY1]=derevo(nt,PXY,DATA.XY,22);

[CrDATA]=CrackProp(DATA,PR.dl);
[nt,PXY]=Tube_perc(PR,CrDATA,DATA.XY,0.9,WXY);

[gt,GS]=Tresh_Gor(1,DATA.XY,PR.Nl);
%nt(:)={[]};
%plot(DATA.XY(nt{:},1),DATA.XY(nt{:},2),'*')
[C,A2C,dVc,pc,DATA.WonV]=Conek(DATA.XY,nt,PR.Nl,CrDATA,DATA.Won,PR.dh,PR.Kc);


%[nt2,PXY2]=derevo(nt,PXY,DATA.XY,23);

gt(:)={[]};
[G,A2G,dVg,pg,DATA.WonG]=Gorizont(DATA.XY,GS,gt,WXY);%Conek(DATA.XY,nt2,PR.Nl,CrDATA,DATA.Won,PR.dh,PR.Kc);%

parfor i=1:8
   SD(i)=Par_Sim_MKT1(PR,DATA,GYData,i,PXY,WXY,Z,C,A2C,G,A2G,dVc,dVg,CrDATA);
end;
end


function CD=Par_Sim_MKT1(PR,DATA,GYData,i,PXY,WXY,Z,C,A2C,G,A2G,dVc,dVg,CrDATA)

[WData]=Well_DATA_P(WXY,Z,PR.Ta,1);
[nt,PXY]=Tube_perc(PR,CrDATA,DATA.XY,1.1,WXY);
[C,A2C,dVc,pc,DATA.WonV]=Conek(DATA.XY,nt,PR.Nl,CrDATA,DATA.Won,PR.dh,PR.Kc);
[XY,KX,Z,Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS]=SimT_MKT(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,1);

SD(1,1)={Q};
SD(2,1)={Sw};
SD(3,1)={Pi};
SD(4,1)={XY};
%SD(5,1)={WXY};
SD(6,1)={MCp};
SD(7,1)={p};
SD(8,1)={PXY};
%SD(9,1)={gt};
SD(10,1)={Z};
SD(11,1)={Ti};
SD(12,1)={Pw};
SD(13,1)={PpW};

CD={SD};

%VZL(XY,K,WXY,Z,Pi,Sw,Cp,Prop.Nl,p,Q,PXY,gt);
we=Q(:,2,:);
i
sum(we(:))
end