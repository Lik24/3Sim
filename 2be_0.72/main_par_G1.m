function CD=main_par_G1(AS)

addpath('Sim_Lib','Tube_Lib','Gor_crack','Sparse_GPU','CrGeom','Termal_lib','GeoMeh_Lib',...
    'Well_lib','Crack_gen','Problems','Poly_lib','SS_lib','Diff_lib','Viz_lib','DATA_In','Adap_lib');
PR=Gl_PRM;%imp_glb_prm;%
[KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z,XY_GY,XY_GY_new,GY_subl]=Sintetic_Real(PR.Ns,PR.Nl);
%dK=rand(size(KX))-1;
% load('dk.mat','dK')
% %save('dk.mat','dK')
% KX=KX+dK.*max(KX(:));
% KY=KY+dK.*max(KY(:));
% KZ=KZ+dK.*max(KZ(:));
SD=load('WHor.mat');
for i=1:9
    WXY=SD.WHO100{i};
    [WData,Won3,~]=Well_DATA(WXY,Z,PR.Ta,PR.Nl,PR.drob);
    [nt,PXY,gXY,PR.dl,tXY,XY_GY1,Won3,WData]=kvad_crack_fun1(XY_GY,PR.Nl,WData,PR.drob,Won3,SD.g_cr_90);
    [DATA]=GridProp(KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY(:,:),H,Z,tXY,PR.Nl,WData.WXY,XY_GY,GY_subl,Won3);
    [GYData,DATA.gKX,DATA.gSw,B,A2B,dVb,pb]=GY_DATA(0,DATA,XY_GY_new,PR); %0/1 - выкл/вкл. аквифер

    DTA(i).WData=WData;
    DTA(i).Won3=Won3;
    DTA(i).nt=nt;
    DTA(i).PR=PR;
    %DTA(i).tXY=tXY;
    DTA(i).XY_GY=XY_GY1;
    DTA(i).DATA=DATA;
    DTA(i).GYData=GYData;
    DTA(i).B=B;
    DTA(i).A2B=A2B;
    DTA(i).dVb=dVb;
    DTA(i).pb=pb;
end



parfor i=1:9
    [CD(:,i),v0(i)]=run_sim(KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,SD.WHO100{i},H,Z,XY_GY_new,GY_subl,DTA(i),i);
end
end
function [CD,V0]=run_sim(KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z,XY_GY_new,GY_subl,DTA,i)
    WData=DTA.WData;
    Won3=DTA.Won3;
    nt=DTA.nt;
    PR=DTA.PR;
    %tXY=DTA.tXY;
    XY_GY=DTA.XY_GY;
    DATA=DTA.DATA;
    GYData=DTA.GYData;
    B=DTA.B;
    A2B=DTA.A2B;
    dVb=DTA.dVb;
    pb=DTA.pb;
    
    Sw0=DATA.gSw;
    
%nt=elka(0,PR.Nl,DATA.XY,15,5,0,25);  %0/1 - выкл/вкл.; кол-во трещин, длинна, флаг к скважине
[CrDATA]=CrackProp(DATA,PR,nt);
%[nt,PXY]=Tube_perc(PR,CrDATA,DATA.XY,1.1,WXY);
[C,A2C,dVc,pc,DATA.WonV,DATA.Lc,CR_GRUP]=Conek2(DATA.XY,nt,PR.Nl,CrDATA,DATA.Won,WData.r0,DATA.ka);

gt=Tresh_Gor(PR.fC(2),DATA.XY,PR.Nl);  % 0/1 - выкл/вкл. горизонтальные трещ.
[GData]=Conek2G(DATA,gt,PR.Nl,CrDATA,WData);

nd=DPorist(PR.fC(3),DATA.XY,PR.Nl); % 0/1 - выкл/вкл. двойная пористость
[DData,~,DATA.gMp]=Conek2D(DATA,nd,PR.Nl,CrDATA,WData,A2C,GData.A2G,PR.ddol);

[Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,Uf,dt1,dV0,DATA.ka,dtz,DSw]=SimT_MKT(PR,C,A2C,GData,B,A2B,DData,dVc,dVb,DATA,WData,GYData,1,CR_GRUP);
i
VZL(DATA,WData.WXY,Pi,Sw(:,end),Ti,MCp,PR.Nl,p,Q,SwC,CR_GRUP,pc,nt,XY_GY,Uf(:,end),pb,GYData,XY_GY_new,dtz,Won3);
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
v0(i)=full(V0);
CD(1)={[Qo,Ql,Qz,sQo,sQl]};
CD(2)={Pw'};
CD(3)={PpW'};
CD(4)={Sw(:,end)};
CD(5)={p };

end