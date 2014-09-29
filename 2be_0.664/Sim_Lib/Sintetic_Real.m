function [KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z,XYc]=Sintetic_Real(Nw,Nl)
SD=load('Buzachi18x19x3.mat');
SD1=load('DATA_WXY.mat');
XYc=[SD.Xc(:),-SD.Yc(:)];
WXY=SD1.WXY;
mXY=min([XYc;WXY]);

XYc(:,1)=XYc(:,1)-mXY(1);
XYc(:,2)=XYc(:,2)-mXY(2);
WXY(:,1)=WXY(:,1)-mXY(1);
WXY(:,2)=WXY(:,2)-mXY(2);


% X=[0,1000,500,0,1000]/4;
% Y=[0,0,500,1000,1000]/4;
% 
% wXY(:,1)=X(1:Nw);  
% wXY(:,2)=Y(1:Nw);
% WXY=wXY;

Mp=SD.Mp;
KX=SD.Kx/1000*8.64;
KY=SD.Ky/1000*8.64;
KZ=SD.Kz/1000*8.64;
Sor=SD.Sor;
Swr=SD.Swr;
Sw=SD.Sw0;

H=SD.DZ.*SD.NTG;
NTG=ones(size(Mp));
n1=size(Mp,1)*size(Mp,2);
n2=size(Mp,3);

Mp=reshape(Mp,n1,n2);
KX=reshape(KX,n1,n2);
KY=reshape(KY,n1,n2);
KZ=reshape(KZ,n1,n2);
Sw=reshape(Sw,n1,n2);
H=reshape(H,n1,n2);
NTG=reshape(NTG,n1,n2);
%KX(:,:)=1;
%Sw(:,2)=1;
Cp=0*ones(size(KX));
T=40*ones(size(KX));
P=100*ones(size(KX))+0.001*rand(size(KX));
z=1:Nl;
Z=repmat(z*20,n1,1);

