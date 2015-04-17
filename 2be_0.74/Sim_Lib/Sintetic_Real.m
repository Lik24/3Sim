function [KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z,XY_GY,XY_GY_new,GY_subl,pXY]=Sintetic_Real(Nw,Nl,fl2)

%[GY,WXY,H,Hk,K,Mp,Sw,Z1,Z3,GY_subl]=read_mr_prop;
%[GY,WXY,H,Hk,K,Mp,Sw,Z1,Z3,GY_subl]=read_mr_prop_MF1;

load('toEmile.mat');
CD=DATA;

Won=CD.Won;
WXY=CD.XY(Won,:);
mXY=min(WXY);

pXY=CD.XY;
pXY(:,1)=pXY(:,1)-mXY(1);
pXY(:,2)=pXY(:,2)-mXY(2);

WXY(:,1)=WXY(:,1)-mXY(1);
WXY(:,2)=WXY(:,2)-mXY(2);
XYmax=max(WXY); Lgr=( sum(XYmax.^2)  )^0.5/2;
Mp=CD.gMp;
K=CD.gKX(:,1);
K(K>10000)=10000;
H=CD.gH(:,1);
Hk=CD.gH(:,1).*CD.gNTG(:,1);

Sw=CD.gSw(:,1);

Cp=CD.gCp(:,1);
T=CD.gT(:,1);
P=CD.gP(:,1);
Z=CD.gZ(:,1);

Mp=repmat(Mp,1,Nl);
%dK=0.5*repmat(K,1,Nl).*(0.5-rand(size(K,1),Nl));
load('dK1')
K=repmat(K,1,Nl)+dK*10;%+0.5*repmat(K,1,Nl).*(0.5-rand(size(K,1),Nl));
H=repmat(H,1,Nl);
Hk=repmat(Hk,1,Nl);
Sw=repmat(Sw,1,Nl);
Cp=repmat(Cp,1,Nl);
T=repmat(T,1,Nl);
P=repmat(P,1,Nl);
Z=repmat(Z,1,Nl)+repmat(1:10:10*Nl,size(Z,1),1);

k=boundary(WXY(:,1),WXY(:,2));
FD=load('GY300');
GY=FD.GY_xy2;%WXY(k,:);
% GY_subl=[[1,-1];GY];
% n6=size(GY_subl,1);
% GY_subl=repmat(GY_subl,Nl,1);
% GY_subl(1:n6:end,1)=1:Nl;

GY1=GY(GY(:,1)<1500,:);
GY3=GY(GY(:,1)>800,:);
GY_subl=[[1,-1];GY1;[2,-1];GY;[3,-1];GY3];

XY_GY=[GY(:,1),GY(:,2)];

SD.Mp=Mp;
SD.K=K;
SD.Sw=Sw;
SD.Z3=Z;
SD.H=H;
SD.Hk=Hk;


pR=300;
min_gy=min(XY_GY);
max_gy=max(XY_GY);
dx=max_gy(1)-min_gy(1);
dy=max_gy(2)-min_gy(2);
xyc=mean(XY_GY(1:end-1,:));
xyc(2)=xyc(2);
if dx>dy
  xc(1)=xyc(1)-0.25*dx;
  xc(2)=xyc(1)+0.25*dx;
  for i=1:size(XY_GY,1)
    if XY_GY(i,1)<xyc(1)
        [t,R]=cart2pol(XY_GY(i,1)-xc(1),XY_GY(i,2)-xyc(2));
        [XY_GY_new(i,1),XY_GY_new(i,2)]=pol2cart(t,R+pR);
        XY_GY_new(i,1)=XY_GY_new(i,1)+xc(1);
        XY_GY_new(i,2)=XY_GY_new(i,2)+xyc(2);
    else
        [t,R]=cart2pol(XY_GY(i,1)-xc(2),XY_GY(i,2)-xyc(2));
        [XY_GY_new(i,1),XY_GY_new(i,2)]=pol2cart(t,R+pR);
        XY_GY_new(i,1)=XY_GY_new(i,1)+xc(2);
        XY_GY_new(i,2)=XY_GY_new(i,2)+xyc(2);
    end
  end;
else
  yc(1)=xyc(2)-0.25*dy;
  yc(2)=xyc(2)+0.25*dy; 
    for i=1:size(XY_GY,1)
    if XY_GY(i,1)<xyc(1)
        [t,R]=cart2pol(XY_GY(i,1)-xyc(1),XY_GY(i,2)-yc(1));
        [XY_GY_new(i,1),XY_GY_new(i,2)]=pol2cart(t,R+pR);
        XY_GY_new(i,1)=XY_GY_new(i,1)+xyc(1);
        XY_GY_new(i,2)=XY_GY_new(i,2)+yc(1);
    else
        [t,R]=cart2pol(XY_GY(i,1)-xyc(1),XY_GY(i,2)-yc(2));
        [XY_GY_new(i,1),XY_GY_new(i,2)]=pol2cart(t,R+pR);
        XY_GY_new(i,1)=XY_GY_new(i,1)+xyc(1);
        XY_GY_new(i,2)=XY_GY_new(i,2)+yc(2);
    end
  end;
end


% [t,R]=cart2pol(XY_GY(:,1)-xyc(1),XY_GY(:,2)-xyc(2));
% [XY_GY_new(:,1),XY_GY_new(:,2)]=pol2cart(t,R+20);
% XY_GY_new(:,1)=XY_GY_new(:,1)+xyc(1);
% XY_GY_new(:,2)=XY_GY_new(:,2)+xyc(2);

%WXY=SD1.WXY;
% WXY2=SD1.WXY2;
% 
% WXY=[WXY;WXY2];

%  WXY4=SD1.WXY4;
%  WXY=[WXY;WXY4];

mXY=min([XY_GY;WXY]);

XY_GY(:,1)=XY_GY(:,1)-mXY(1);
XY_GY(:,2)=XY_GY(:,2)-mXY(2);
XY_GY_new(:,1)=XY_GY_new(:,1)-mXY(1);
XY_GY_new(:,2)=XY_GY_new(:,2)-mXY(2);
WXY(:,1)=WXY(:,1)-mXY(1);
WXY(:,2)=WXY(:,2)-mXY(2);

GY_subl1=GY_subl;
GY_subl(:,1)=GY_subl(:,1)-mXY(1);
GY_subl(:,2)=GY_subl(:,2)-mXY(2);
GY_subl(GY_subl1(:,2)==-1,:)=GY_subl1(GY_subl1(:,2)==-1,:);

% X=[0,1000,500,0,1000]/4;
% Y=[0,0,500,1000,1000]/4;
% 
% wXY(:,1)=X(1:Nw);  
% wXY(:,2)=Y(1:Nw);
% WXY=wXY;

Mp=SD.Mp;
KX=SD.K/1000*8.64;
KY=SD.K/1000*8.64;
KZ=SD.K/1000*8.64*0.1*0;

Sw=SD.Sw;
Z=SD.Z3;

H=SD.H;
NTG=SD.Hk./H;
NTG(H==0)=0;
n1=size(Mp,1);
n2=size(Mp,2);

Mp=reshape(Mp,n1,n2);
KX=reshape(KX,n1,n2);
KY=reshape(KY,n1,n2);
KZ=reshape(KZ,n1,n2);
Sw=reshape(Sw,n1,n2);
H=reshape(H,n1,n2);
NTG=reshape(NTG,n1,n2);
%KX(:,:)=1;
%Sw(:,2)=1;
% for l=1:n2
%  Sw(:,l)=mean(Sw(:,l));
% end
Cp=0*ones(size(KX));
T=40*ones(size(KX));
%P=100*ones(size(KX))+0.001*rand(size(KX));
P=9.81*1000*Z*1e-5;
z=1:Nl;
