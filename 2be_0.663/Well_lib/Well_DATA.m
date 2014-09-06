function [WData]=Well_DATA(WXY,Z,Ta)
r0=0.084;

n=size(WXY);
n1=size(Z);
%load('DATA_WXY','a0','a1','a2')

%a=rand(n(1),1);
  Uf=-ones(4,1);
  Uf(1)=1;
  %Uf(2:3)=0;
 % Uf=a0;

%Uf=uf2;

Pw(Uf==1,1)=50;
Pw(Uf==-1,1)=150;

Qz(Uf(1:n(1))==1,1)=0;
Qz(Uf(1:n(1))==-1,1)=-5*0;

Uf=repmat(Uf,n1(2),1);

CpW=zeros(n(1),1);
CpW(Uf==1,1)=0;
CpW(Uf==-1,1)=0;

TeW=40+zeros(n(1),1);
TeW(Uf(1:n(1))==-1,1)=100;

%Uf([1,3:4,6])=0;
%Uf([10:27])=0;

WData.Pw=repmat(Pw,1,Ta);
WData.Uf=repmat(Uf,1,Ta);
WData.CpW=repmat(CpW,1,Ta);
WData.Qz=repmat(Qz,1,Ta);
WData.TeW=repmat(TeW,1,Ta);
WData.Doly=ones(n(1),1);
WData.r0=r0;
WData.WXY=WXY;

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