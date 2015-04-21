function SD=main_par(fl,DATA,GYData,PR,WData,DIN)

addpath('Sim_Lib','Tube_Lib','Gor_crack','Sparse_GPU','CrGeom','Termal_lib','GeoMeh_Lib',...
    'Well_lib','Crack_gen','Problems','Poly_lib','SS_lib','Diff_lib','Viz_lib','DATA_In','Adap_lib');

parfor i=1:4
    
   SD(i)=par_run(DATA(i),PR,fl,GYData,WData,DIN(:,i));
end
end

function SD=par_run(DATA,PR,fl,GYData,WData,DIN)

B=DIN{1};
A2B=DIN{2};
dVb=DIN{3};
pb=DIN{4};
nt=DIN{5};
% gXY=DIN(6);
% tXY=DIN(7);
% XY_GY=DIN{8};
% Won3=DIN{9};

Sw0=DATA.gSw;
if fl(2)==0
nt={zeros(2,0)};
end
%nt=elka(1,PR.Nl,DATA.XY,3,30,0,25);  %0/1 - выкл/вкл.; кол-во трещин, длинна, флаг к скважине

[CrDATA]=CrackProp(DATA,PR,nt);
%[nt,PXY]=Tube_perc(PR,CrDATA,DATA.XY,1.1,WXY);
[C,A2C,dVc,pc,DATA.WonV,DATA.Lc,CR_GRUP]=Conek2(DATA.XY,DATA.gZ,nt,PR.Nl,CrDATA,DATA.Won,WData.r0,DATA.ka);

gt=Tresh_Gor(PR.fC(2),DATA.XY,PR.Nl);  % 0/1 - выкл/вкл. горизонтальные трещ.
[GData]=Conek2G(DATA,gt,PR.Nl,CrDATA,WData);

nd=DPorist(PR.fC(3),DATA.XY,PR.Nl); % 0/1 - выкл/вкл. двойная пористость
[DData,~,DATA.gMp]=Conek2D(DATA,nd,PR.Nl,CrDATA,WData,A2C,GData.A2G,PR.ddol);

[Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,Uf,dt1,dV0,DATA.ka,dtz,DSw]=SimT_MKT(PR,C,A2C,GData,B,A2B,DData,dVc,dVb,DATA,WData,GYData,1,CR_GRUP);

%VZL(DATA,WData.WXY,Pi,Sw(:,end),Ti,MCp,PR.Nl,p,Q,SwC,CR_GRUP,pc,nt,XY_GY,Uf(:,end),pb,GYData,XY_GY_new,dtz,Won3,A2C);
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