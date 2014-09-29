function [KX,KY,KZ,Mp,P,Sw,Cp,T,NTG,WXY,H,Z]=Sintetic(Nw,Nl)
%WXY=rand(Nw,2)*1000;
%[X,Y]=meshgrid(100:400:1000,100:300:1000);
X=[0,1000,0,1000]/4;
Y=[0,0,1000,1000]/4;

% X=[0,500,1000,0,500,1000,0,500,1000];
% Y=[0,0,0,500,500,500,1000,1000,1000];

WXY(:,1)=X(1:Nw);  
WXY(:,2)=Y(1:Nw);

% v=convhull(WXY);
% S=polyarea(WXY(v,1),WXY(v,2));
KX=rand(Nw,Nl);
KY=rand(Nw,Nl);
KZ=rand(Nw,Nl);
%KX
KX(:,:)=1.0*8.64;%mean(K);
KX(:,2)=1.0*8.64;
%KY
KY(:,:)=1.0*8.64;%mean(K);
KY(:,2)=1.0*8.64;
%KZ
KZ(:,:)=0.01*8.64;%mean(K);
KZ(:,2)=0.01*8.64;

H=10*ones(size(KX));
%H(:,2)=7.5*ones(size(K,1),1);
Mp=0.2*ones(size(KX));
NTG=1*ones(size(KX));

Sw=0*ones(size(KX));
%Sw(:,2)=1;
Cp=0.0*ones(size(KX));
T=40*ones(size(KX));
P=100*ones(size(KX));
z=1:Nl;
Z=repmat(z*20,Nw,1);

