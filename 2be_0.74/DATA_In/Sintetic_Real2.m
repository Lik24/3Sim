function [KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z,Nw,Lgr]=Sintetic_Real2
load('toEmile.mat');
SD=DATA;
Won=SD.Won;
%SD1=load('WCoord.mat');
%XYc=[SD.Xc(:),-SD.Yc(:)];
WXY=SD.XY(Won,:);
mXY=min(WXY);

%XYc(:,1)=XYc(:,1)-mXY(1);
%XYc(:,2)=XYc(:,2)-mXY(2);
WXY(:,1)=WXY(:,1)-mXY(1);
WXY(:,2)=WXY(:,2)-mXY(2);
XYmax=max(WXY); Lgr=( sum(XYmax.^2)  )^0.5/2;
Mp=SD.gMp;
KX=SD.gKX(Won,1);
KY=SD.gKY(Won,1);
KZ=SD.gKZ(Won,1);
H=SD.gH(Won,1);
Nw=max(Won);
Sw=SD.gSw(Won,1);
NTG=SD.gNTG(Won,1);
Cp=SD.gCp(Won,1);
T=SD.gT(Won,1);
P=SD.gP(Won,1);
Z=SD.gZ(Won,1);

% X=[0,1000,500,0,1000]/4;
% Y=[0,0,500,1000,1000]/4;
% 
% wXY(:,1)=X(1:Nw);  
% wXY(:,2)=Y(1:Nw);
% WXY=wXY;


% ImgMp=zeros(size(Mp0)); 
% ImgMp(find(Mp0~=0))=1;
% Mp=sum(Mp0,2)./sum(ImgMp,2);


% ImgK=zeros(size(K0));
% ImgK(find(K0~=0))=1;
% KX=sum(K0,2)./sum(ImgK,2);

%KY=SD.Ky/1000*8.64;
%KZ=SD.Kz/1000*8.64;
% Sor=SD.Sor;
% Swr=SD.Swr;



% ImgH=zeros(size(H0));
% ImgH(find(H0~=0))=1;
% H=sum(H0,2)./sum(ImgH,2);

% kk=(sum(ImgH,2).*sum(ImgK,2).*sum(ImgMp,2)>0);

% Mp=Mp(kk,1);
% KX=KX(kk,1);
% H=H(kk,1);
% 
% KY=0*ones(size(KX))*8.64;
% KZ=0*ones(size(KX))*8.64;

% WXY=WXY(kk,:);

% n1=size(Mp,1)*size(Mp,2);
% n2=size(Mp,3);
% 
% Mp=reshape(Mp,n1,n2);
% KX=reshape(KX,n1,n2);
% KY=reshape(KY,n1,n2);
% KZ=reshape(KZ,n1,n2);
% Sw=reshape(Sw,n1,n2);
% H=reshape(H,n1,n2);
% NTG=reshape(NTG,n1,n2);
%KX(:,:)=1;
%Sw(:,2)=1;

end

