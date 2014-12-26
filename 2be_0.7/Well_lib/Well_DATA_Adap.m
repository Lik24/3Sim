function [WData,Ppl_d,Pw_d,Pw_z]=Well_DATA_Adap(WXY,Z,Ta)
r0=0.084;

n=size(WXY);
n1=size(Z);
load('history_deb_day','a0','dob_lik2','zak_lik2','dob_Ppl2','dob_Pw2','zak_Pw2');
Ppl_d=dob_Ppl2;
Pw_d=dob_Pw2;
Pw_z=zak_Pw2;

Pw_d(isnan(Pw_d)==1)=0;
Pw_z(isnan(Pw_z)==1)=0;

%a=rand(n(1),1);
  Uf=ones(33,1);
 % Uf(1:3:end)=-1;
  Uf=a0(:,end);

%Uf=uf2;

Pw(Uf==1,1)=34;
Pw(Uf==-1,1)=150;

Qz(Uf(1:n(1))==1,1)=0;
Qz(Uf(1:n(1))==-1,1)=-15;
%Qz(:,1)=dob_lik2-zak_lik2;

Uf=repmat(Uf,n1(2),1);

CpW=zeros(n(1),1);
CpW(Uf==1,1)=0;
CpW(Uf==-1,1)=0;

TeW=40+zeros(n(1),1);
TeW(Uf(1:n(1))==-1,1)=100;

Pw_Q_C_bnd=[-200,700;   % Ограничения на Забойные давления
            0,1000;     % Ограничения на дебит доб. скв.
            -1000,0;    % Ограничения на дебит наг. скв.
            0.98,0.98;  % Предельная обв.
            0,0];       % Предельный дебит нефти.
  
WData.Pw=repmat(Pw,1,Ta);
WData.Uf=repmat(Uf,1,Ta);
WData.CpW=repmat(CpW,1,Ta);
WData.Qz=repmat(Qz,1,Ta);
WData.TeW=repmat(TeW,1,Ta);
WData.Doly=ones(n(1),n1(2));
WData.r0=r0;
WData.WXY=WXY;
WData.PwQC_bnd=Pw_Q_C_bnd;

 WData.Uf=a0;
 WData.Qz=dob_lik2-zak_lik2;
 WQ=WData.Qz;
 for i=1:size(WData.Qz,1)
 WData.Qz(i,:)=smooth(WData.Qz(i,:));
 WData.Qz(i,WQ(i,:)==0)=0;
 end;

