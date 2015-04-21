function [WData]=Well_DATA_mP(WXY,Z,Ta,i)
r0=0.084;

n=size(WXY);
n1=size(Z);
load('history_deb_day','a0','dob_lik2','zak_lik2');
AS=load('WCoord');
SD=load('ADP_ReZ2');
a0_h=a0;
a0=[a0_h(:,end);AS.a04];
ED=load('a0_in');
a0=ED.a0_sol2(:,i);
a0_old=ED.a0_sol(:,1);

%a=rand(n(1),1);
  Uf=ones(33,1);
  %Uf(1:3:end)=-1;
  Uf=a0(:,end);
%  Uf(Uf<0)=0;
%Uf=uf2;

%Pw(Uf==1,1)=2;
%Pw(Uf==-1,1)=150;
Pw_in=mean(SD.Pw0(a0_h(:,end)==-1));
Pw_dob=mean(SD.Pw0(a0_h(:,end)==1));

Pw=zeros(size(a0));
Pw(a0==-1,1)=Pw_in;
Pw(a0==1,1)=Pw_dob;

Pw_old=SD.Pw0(:,end);
Pw_old=Pw_old.*(a0(1:33)==a0_old(1:33));
a0_new=a0(1:33);
Pw_old(Pw_old==0)=Pw_in.*(a0_new(Pw_old==0)==-1)+Pw_dob.*(a0_new(Pw_old==0)==1);
Pw(1:33,1)=Pw_old;

Qz=zeros(size(Uf(1:n(1))));
Qz(Uf(1:n(1))==1,1)=0;
Qz(Uf(1:n(1))==-1,1)=-100*(i~=2)-75*(i==2);
%Qz(:,1)=dob_lik2-zak_lik2;

Uf=repmat(Uf,n1(2),1);

CpW=zeros(n(1),1);
CpW(Uf==1,1)=0;
CpW(Uf==-1,1)=0;

TeW=40+zeros(n(1),1);
TeW(Uf(1:n(1))==-1,1)=100;

Pw_Q_C_bnd=[1,70;   % Ограничения на Забойные давления
            0,100;     % Ограничения на дебит доб. скв.
            -100,0;    % Ограничения на дебит наг. скв.
            0.98,0.98;  % Предельная обв.
            1,1];       % Предельный дебит нефти.

WData.Pw=repmat(Pw,1,Ta);
WData.Uf=repmat(Uf,1,Ta);
WData.CpW=repmat(CpW,1,Ta);
WData.Qz=repmat(Qz,1,Ta);
WData.TeW=repmat(TeW,1,Ta);
WData.Doly=ones(n(1),n1(2));
WData.r0=r0;
WData.WXY=WXY;
WData.PwQC_bnd=Pw_Q_C_bnd;

%WData.Uf=a0;
%WData.Qz=dob_lik2-zak_lik2;
% for i=1:Ta
%  if mod(i,8*280)<8*140
%    WData.Uf(:,i)=uf1(:,i);  
%  else
%    WData.Uf(:,i)=uf2(:,i);  
%  end;
%   uf=WData.Uf(1:n(1),i);
%  WData.Qz(uf==-1,i)=-125; 
% end;

%WData.CpW(:,ceil(end/3*2):end)=0.01+WData.CpW(:,ceil(end/3*2):end);
% plot(sum(WData.Uf,1))
% %plot(rt)
% fggh