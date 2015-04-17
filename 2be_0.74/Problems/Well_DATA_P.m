function [WData]=Well_DATA_P(WXY,Z,Ta,i)
r0=0.084;

n=size(WXY);
n1=size(Z);
load('DATA_WXY','a0','a1','a2')

%a=rand(n(1),1);
Uf=a0;
uf1=Uf;
uf2=Uf;

uf1(a1)=0;
uf2(a2)=0;

if i==2
Uf=uf1;
elseif i==3
Uf=uf2;
end;


uf1=repmat(uf1,n1(2),1);
uf2=repmat(uf2,n1(2),1);
uf1=repmat(uf1,1,Ta);
uf2=repmat(uf2,1,Ta);


Pw(Uf==1,1)=50;
Pw(Uf==-1,1)=150;   

Qz(Uf(1:n(1))==1,1)=0;

if i==1
    Qz(Uf(1:n(1))==-1,1)=-125;
elseif i==2
    Qz(Uf(1:n(1))==-1,1)=-125*2;
elseif i==3
    Qz(Uf(1:n(1))==-1,1)=-125*2;
end;

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
WData.Doly=ones(n(1)*3,1);
WData.r0=r0;
WData.WXY=WXY;

dt=[1,1,1,10,60,365,2*365,4*365];
if i>3
    dT=dt(i);
    for j=365*10:Ta
        if mod(j,dT)<dT/2
            WData.Uf(:,j)=uf1(:,j);
        else
            WData.Uf(:,j)=uf2(:,j);
        end;
        uf=WData.Uf(1:n(1),j);
        WData.Qz(uf==-1,j)=-125*2;
        WData.Qz(uf~=-1,j)=0;
    end;
end;

%WData.CpW(:,ceil(end/3*2):end)=0.01+WData.CpW(:,ceil(end/3*2):end);
% plot(sum(WData.Uf,1))
% %plot(rt)
% fggh