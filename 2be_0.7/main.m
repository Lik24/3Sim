clearvars

addpath('Sim_Lib','Tube_Lib','Gor_crack','Sparse_GPU','CrGeom','Termal_lib','GeoMeh_Lib',...
    'Well_lib','Crack_gen','Problems','Poly_lib','SS_lib','Diff_lib','Viz_lib','DATA_In','Adap_lib','2exe');
PR=imp_glb_prm;%Gl_PRM;

[KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z,XY_GY,XY_GY_new,GY_subl]=Sintetic_Real(PR.Ns,PR.Nl);
%Sw(:)=0.3;
%KX(:)=mean(KX(:));
XY_GYs=XY_GY;
%[KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z]=Sintetic(PR.Ns,PR.Nl);
[WData]=Well_DATA(WXY,Z,PR.Ta,PR.Nl);
%[WData,Ppwf,Pw_d,Pw_z]=Well_DATA_Adap(WXY,Z,PR.Ta);

[nt,PXY,gXY,PR.dl,tXY,XY_GY]=kvad_crack_fun(XY_GY,PR.Nl,WXY);
[DATA]=GridProp(KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY(:,:),H,Z,tXY,PR.Nl,WXY,XY_GY,GY_subl);
Sw0=DATA.gSw;
[GYData,DATA.gKX,DATA.gSw,B,A2B,dVb,pb]=GY_DATA(0,DATA,XY_GY_new,PR); %0/1 - выкл/вкл. аквифер

%[WData.Doly,DATA.gKX,GYData.GY_Kxy]=Load_adp_prm2(WData.Doly,DATA.gKX,GYData.GY_Kxy);
%[WData.Doly,DATA,GYData]=Load_adp_prm(DATA,GYData,tXY);
%[nt1,PXY]=derevo(nt,DATA.XY,22);

%nt=elka(0,PR.Nl,DATA.XY,6,10,0,25);  %0/1 - выкл/вкл.; кол-во трещин, длинна, флаг к скважине
[CrDATA]=CrackProp(DATA,PR,nt);
%[nt,PXY]=Tube_perc(PR,CrDATA,DATA.XY,1.1,WXY);
[C,A2C,dVc,pc,DATA.WonV,DATA.Lc,CR_GRUP]=Conek2(DATA.XY,nt,PR.Nl,CrDATA,DATA.Won,WData.r0,DATA.ka);

gt=Tresh_Gor(0,DATA.XY,PR.Nl);  % 0/1 - выкл/вкл. горизонтальные трещ.
[GData]=Conek2G(DATA,gt,PR.Nl,CrDATA,WData);

nd=DPorist(1,DATA.XY,PR.Nl); % 0/1 - выкл/вкл. двойная пористость
[DData,~,DATA.gMp]=Conek2D(DATA,nd,PR.Nl,CrDATA,WData,A2C,GData.A2G);

[Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,Uf,dt1,dV0,DATA.ka,dtz,DSw]=SimT_MKT(PR,C,A2C,GData,B,A2B,DData,dVc,dVb,DATA,WData,GYData,1,CR_GRUP);

VZL(DATA,WXY,Pi,Sw(:,end),Ti,MCp,PR.Nl,p,Q,SwC,CR_GRUP,pc,nt,XY_GY,Uf(:,end),pb,GYData,XY_GY_new,dtz);
%VZL_VORONOI(DATA,Pi(:,end),p,WXY,WData.Uf(:,end))

  Qo(:,1)=sum(Q(:,3,:));
  Ql(:,1)=sum(Q(:,2,:));
  Qz(:,1)=sum(Q(:,1,:));

 sQo=cumsum(Qo,1);
 sQl=cumsum(Ql,1);

c=1-Qo./Ql;
%plot(c)
dV1([p,size(p,2)+DData.pd])=dV0([1:size(p,2),size(p,2)+size(pc,2)+1:size(pc,2)+size(p,2)+size(DData.pd,2)]);
Sw0=Sw0(DATA.ka==1); Sw0=[Sw0;Sw0(1:size(DData.D,1))];
V0=sum(dV1'.*(1-Sw0));
sQo(end,:)/V0

TBL=table(single(Ql),single(Qo),single(Qz),single(c),'VariableNames',{'Ql','Qo','Qz','c'});
writetable(TBL,'OutQ.txt','Delimiter','tab')