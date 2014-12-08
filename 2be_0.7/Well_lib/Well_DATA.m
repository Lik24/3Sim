function [WData]=Well_DATA(WXY,Z,Ta)
r0=0.084;

n=size(WXY);
n1=size(Z);
%n1(2)=2;
%load('history_deb_day_ura','a0','dob_lik2','zak_lik2');

  %a=rand(n(1),1);
  Uf=ones(4,1);
  Uf(1:10:end)=-1;
  %Uf(1:3:end)=-1;
  %Uf=a0(:,end);
%  Uf(Uf<0)=0;
%Uf=uf2;

Pw(Uf==1,1)=2;
Pw(Uf==-1,1)=170;
% Pw_in=mean(SD.Pw0(a0_h(:,end)==-1));
% Pw_dob=mean(SD.Pw0(a0_h(:,end)==1));


Qz(Uf(1:n(1))==1,1)=0;
Qz(Uf(1:n(1))==-1,1)=-150;
%Qz(:,1)=dob_lik2-zak_lik2;

Uf=repmat(Uf,n1(2),1);

CpW=zeros(n(1),1);
CpW(Uf==1,1)=0;
CpW(Uf==-1,1)=0;

TeW=40+zeros(n(1),1);
TeW(Uf(1:n(1))==-1,1)=100;

Pw_Q_C_bnd=[1,70;   % Ограничения на Забойные давления
            0,200;     % Ограничения на дебит доб. скв.
            -200,0;    % Ограничения на дебит наг. скв.
            0.98,0.98;  % Предельная обв.
            1,1];       % Предельный дебит нефти.

WData.Pw=repmat(Pw,1,Ta);
WData.Uf=repmat(Uf,1,Ta);
WData.CpW=repmat(CpW,1,Ta);
WData.Qz=repmat(Qz,1,Ta);
WData.TeW=repmat(TeW,1,Ta);
WData.Doly=ones(n(1),n1(2));  %Коэф связи скважины с пластом
WData.SDoly=ones(n(1),1);     %Коэф связи скважины с меторождением
WData.r0=r0;
WData.WXY=WXY;
WData.PwQC_bnd=Pw_Q_C_bnd;
WData.Uf=repmat(Uf,1,Ta);
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
% fgg