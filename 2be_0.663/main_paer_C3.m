function [SD,Ql,Qo,Qz,c,sQo]=main_paer_C3

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
%[nt1,PXY1]=derevo(nt,DATA.XY,22);

[CrDATA]=CrackProp(DATA,PR.dl);

%[nt,PXY]=Tube_perc(PR,CrDATA,DATA.XY,1.1,WXY);

[gt,GS]=Tresh_Gor(1,DATA.XY,PR.Nl);

% FG=load('derevoNT2.mat','nt1','nt2','PXY1','PXY2');
% nt1=FG.nt1;
% nt2=FG.nt2;
% PXY1=FG.PXY1;
% PXY2=FG.PXY2;
nt(:)={[]};
[C,A2C,dVc,pc,DATA.WonV]=Conek(DATA.XY,nt,PR.Nl,CrDATA,DATA.Won,PR.dh,PR.Kc,WData.r0);

%[nt2,PXY2]=derevo(nt,DATA.XY,24);

gt(:)={[]};
[G,A2G,dVg,pg,DATA.WonG]=Gorizont(DATA.XY,GS,gt,WXY);%Conek(DATA.XY,nt2,PR.Nl,CrDATA,DATA.Won,PR.dh,PR.Kc,WData.r0);%
%save('derevoNT2.mat','nt1','nt2','PXY1','PXY2')
parfor i=1:8
   SD(i)=Par_Sim_MKT1(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,i,PXY,PXY,Z,WXY);
end;

for i=1:size(SD,2)
  CD=SD{i};
  Q=CD{1};
  Qo(:,i)=sum(Q(:,3,:));
  Ql(:,i)=sum(Q(:,2,:));
  Qz(:,i)=sum(Q(:,1,:));
end;

for i=1:size(Qo,1)
 sQo(i,:)=sum(Qo(1:i,:),1);
 sQl(i,:)=sum(Ql(1:i,:),1);
end;

c=1-Qo./Ql;

end


function CD=Par_Sim_MKT1(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,i,PXY1,PXY2,Z,WXY)
[WData]=Well_DATA_P2(WXY,Z,PR.Ta,i);
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