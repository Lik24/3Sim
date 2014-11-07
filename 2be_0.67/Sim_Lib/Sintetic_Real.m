function [KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z,XY_GY,XY_GY_new]=Sintetic_Real(Nw,Nl)
SD=load('RIGIS_DATA_A1.mat');
SD1=load('WCoord.mat');
XY_GY=[SD1.GY(:,1),SD1.GY(:,2)];

pR=200;
min_gy=min(XY_GY);
max_gy=max(XY_GY);
dx=max_gy(1)-min_gy(1);
dy=max_gy(2)-min_gy(2);
xyc=mean(XY_GY(1:end-1,:));
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

WXY=SD1.WXY;
% WXY2=SD1.WXY2;
% 
% WXY=[WXY;WXY2];

%  WXY4=SD1.WXY4;
%  WXY=[WXY;WXY4];

mXY=min([XY_GY;XY_GY_new;WXY]);

XY_GY(:,1)=XY_GY(:,1)-mXY(1);
XY_GY(:,2)=XY_GY(:,2)-mXY(2);
XY_GY_new(:,1)=XY_GY_new(:,1)-mXY(1);
XY_GY_new(:,2)=XY_GY_new(:,2)-mXY(2);
WXY(:,1)=WXY(:,1)-mXY(1);
WXY(:,2)=WXY(:,2)-mXY(2);


% X=[0,1000,500,0,1000]/4;
% Y=[0,0,500,1000,1000]/4;
% 
% wXY(:,1)=X(1:Nw);  
% wXY(:,2)=Y(1:Nw);
% WXY=wXY;

Mp=SD.Mp;
KX=SD.K/1000*8.64;
KY=SD.K/1000*8.64;
KZ=SD.K/1000*8.64;
Sor=0;%SD.Sor;
Swr=0;%SD.Swr;
Sw=SD.Sw;
Z=SD.Z;

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
