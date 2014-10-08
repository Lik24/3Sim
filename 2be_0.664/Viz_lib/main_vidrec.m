
%load('matlab16.11.mat');
%load('derevoNT9.mat','PXY1')

%SD=SD;
Nl=1;
CD=SD{7};
Q=CD{1,1};
Sw=CD{2,1};
Pi=CD{3,1};
XY=CD{4,1};
WXY=CD{5,1};
Cp=CD{6,1};
p=CD{7,1};
PXY=CD{8,1};
%PXY2=PXY;

gt=CD{9,1};
Z=CD{10,1};
uf=CD{11,1};
WXY=zeros(4,2);
WXY(:,1)=[0,250,0,250]';
WXY(:,2)=[0,0,250,250]';
uf=[1,-1,-1,-1]';
video_save(XY,WXY,Z,Pi(:,1:1:end),Sw(:,1:1:end),Nl,p,PXY,PXY,uf)