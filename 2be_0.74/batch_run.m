fl=[1,1,1,5];

addpath('Sim_Lib','Tube_Lib','Gor_crack','Sparse_GPU','CrGeom','Termal_lib','GeoMeh_Lib',...
    'Well_lib','Crack_gen','Problems','Poly_lib','SS_lib','Diff_lib','Viz_lib','DATA_In','Adap_lib');

PR=Gl_PRM;%imp_glb_prm;%
if fl(4)==9
    PR.Ta=60*365;
end
[KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z,XY_GY,XY_GY_new,GY_subl]=Sintetic_Real(PR.Ns,PR.Nl,fl(3));

if fl(1)==1
sD=load('dK16_10');
KX=sD.dK.*KX;
KY=sD.dK.*KY;
KZ=sD.dK.*KZ;
end

KZ0=KZ;
kz(1)={KZ};
%KZ(:,1:8:end)=KZ(:,1:8:end)*0.001; kz(2)={KZ};
KZ(:,1:4:end)=KZ0(:,1:4:end)*0.01; kz(2)={KZ};
KZ(:,1:2:end)=KZ0(:,1:2:end)*0.01; kz(3)={KZ};
KZ(:,1:1:end)=KZ0(:,1:1:end)*0.01; kz(4)={KZ};

for i=1:4
[WData,Won3,~]=Well_DATA(WXY,Z,PR.Ta,PR.Nl,PR.drob,fl(4));
[nt,gXY,PR.dl,tXY,~,Won3,WData]=kvad_crack_fun(XY_GY,PR.Nl,WData,PR.drob,Won3);
[DATA(i)]=GridProp(KX,KY,kz{i},Mp,P,Sw,Cp,T,NTG,WXY(:,:),H,Z,tXY,PR.Nl,WData.WXY,XY_GY,GY_subl,Won3);
[GYData,DATA(i).gKX,DATA(i).gSw,B,A2B,dVb,pb]=GY_DATA(0,DATA(i),XY_GY_new,PR); %0/1 - выкл/вкл. аквифер
DIN(1,i)={B};
DIN(2,i)={A2B};
DIN(3,i)={dVb};
DIN(4,i)={pb};
DIN(5,i)={nt};
% DIN(6,i)={gXY};
% DIN(7,i)={tXY};
% DIN(8,i)={XY_GY};
% DIN(9,i)={Won3};

end

%main_par(fl,DATA,GYData,PR,WData,DIN)
%j3=batch(@main_par,1,{fl,DATA,GYData,PR,WData,DIN},'Profile','local','Pool',3);
j6=batch(@main_par,1,{fl,DATA,GYData,PR,WData,DIN},'Profile','itpm_ccord','Pool',4);
clear DATA GYData WData KX KY KZ