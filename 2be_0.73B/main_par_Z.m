function SD=main_par_Z

addpath('Sim_Lib','Tube_Lib','Gor_crack','Sparse_GPU','CrGeom','Termal_lib','GeoMeh_Lib',...
    'Well_lib','Crack_gen','Problems','Poly_lib','SS_lib','Diff_lib','Viz_lib','DATA_In','Adap_lib');
PR=Gl_PRM;%imp_glb_prm;%

[KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z,XY_GY,XY_GY_new,GY_subl]=Sintetic_Real(PR.Ns,PR.Nl);
kz(1)={KZ};
KZ(:,8:8:end)=KZ(:,8:8:end)*0.1; kz(2)={KZ};
KZ(:,4:4:end)=KZ(:,4:4:end)*0.1; kz(3)={KZ};
KZ(:,2:2:end)=KZ(:,2:2:end)*0.1; kz(4)={KZ};
KZ(:,1:1:end)=KZ(:,1:1:end)*0.1; kz(5)={KZ};

for i=1:5
 [WData,Won3,~]=Well_DATA(WXY,Z,PR.Ta,PR.Nl,PR.drob);
 [nt,gXY,PR.dl,tXY,XY_GY1,Won31,WData1]=kvad_crack_fun(XY_GY,PR.Nl,WData,PR.drob,Won3);
 [DATA(i)]=GridProp(KX,KY,kz{i},Mp,P,Sw,Cp,T,NTG,WXY(:,:),H,Z,tXY,PR.Nl,WData1.WXY,XY_GY1,GY_subl,Won31);

 [GYData(i),DATA(i).gKX,DATA(i).gSw,B,A2B,dVb,pb]=GY_DATA(0,DATA(i),XY_GY_new,PR); %0/1 - выкл/вкл. аквифер
 DT(i).B=B;
 DT(i).A2B=A2B;
 DT(i).dVb=dVb;
 DT(i).pb=pb;
 
end

parfor i=1:5
    SD(i)=par_run(DATA(i),XY_GY,XY_GY_new,PR,nt,Won3,WData1,GYData(i),DT(i));
end
end

function SD=par_run(DATA,XY_GY,XY_GY_new,PR,nt,Won3,WData,GYData,DT)

 B=DT.B;
 A2B=DT.A2B;
 dVb=DT.dVb;
 pb=DT.pb;
 Sw0=DATA.gSw;
  
[CrDATA]=CrackProp(DATA,PR,nt);
%[nt,PXY]=Tube_perc(PR,CrDATA,DATA.XY,1.1,WXY);
[C,A2C,dVc,pc,DATA.WonV,DATA.Lc,CR_GRUP]=Conek2(DATA.XY,DATA.gZ,nt,PR.Nl,CrDATA,DATA.Won,WData.r0,DATA.ka);

gt=Tresh_Gor(PR.fC(2),DATA.XY,PR.Nl);  % 0/1 - выкл/вкл. горизонтальные трещ.
[GData]=Conek2G(DATA,gt,PR.Nl,CrDATA,WData);

nd=DPorist(PR.fC(3),DATA.XY,PR.Nl); % 0/1 - выкл/вкл. двойная пористость
[DData,~,DATA.gMp]=Conek2D(DATA,nd,PR.Nl,CrDATA,WData,A2C,GData.A2G,PR.ddol);

[Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,Uf,dt1,dV0,DATA.ka,dtz,DSw]=SimT_MKT(PR,C,A2C,GData,B,A2B,DData,dVc,dVb,DATA,WData,GYData,1,CR_GRUP);

VZL(DATA,WData.WXY,Pi,Sw(:,end),Ti,MCp,PR.Nl,p,Q,SwC,CR_GRUP,pc,nt,XY_GY,Uf(:,end),pb,GYData,XY_GY_new,dtz,Won3,A2C);
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
sQo(end,:)/V0;
SD={[c,sQo,sQl,sQo/(500*500*40*0.3*0.75),Qo,Ql]};
end