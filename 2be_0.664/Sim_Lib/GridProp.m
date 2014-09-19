function [DATA]=GridProp(wKx,wKy,wKz,wMp,wP,wSw,wCp,wT,wNTG,CXY,wH,wZ,gXY,nl,WXY)
[nw]=size(wKx,1);
A=zeros(nw,nl*6);
XY=[WXY;gXY];

XY=NODuble(XY);
% DT = delaunayTriangulation(XY(:,1),XY(:,2));
% triplot(DT)
% hold on
ncg=size(XY,1);

BndXY=zeros(ncg,1);
ao=convhull(XY(:,1),XY(:,2));
BndXY(ao)=1;
BndXY=repmat(BndXY,nl,1);

gKX=zeros(size(XY,1),nl);
gKY=zeros(size(XY,1),nl);
gKZ=zeros(size(XY,1),nl);
gMp=zeros(size(XY,1),nl);
gP=zeros(size(XY,1),nl);
gSw=zeros(size(XY,1),nl);
gCp=zeros(size(XY,1),nl);
gH=zeros(size(XY,1),nl);
gZ=zeros(size(XY,1),nl);
gT=zeros(size(XY,1),nl);
gNTG=zeros(size(XY,1),nl);

mo=11;
for i=1:nl
    A(:,1+mo*(i-1))=wKx(:,i);
    A(:,2+mo*(i-1))=wMp(:,i);
    A(:,3+mo*(i-1))=wP(:,i);
    A(:,4+mo*(i-1))=wSw(:,i);
    A(:,5+mo*(i-1))=wH(:,i);
    A(:,6+mo*(i-1))=wZ(:,i);
    A(:,7+mo*(i-1))=wCp(:,i);
    A(:,8+mo*(i-1))=wT(:,i);
    A(:,9+mo*(i-1))=wNTG(:,i);
    A(:,10+mo*(i-1))=wKy(:,i);
    A(:,11+mo*(i-1))=wKz(:,i);
end;

parfor i=1:mo*nl
Aa(:,i)=razmaz(A(:,i),CXY,XY);
end

for i=1:nl
    gKX(:,i)=Aa(:,1+mo*(i-1));
    gMp(:,i)=Aa(:,2+mo*(i-1));
    gP(:,i)=Aa(:,3+mo*(i-1));
    gSw(:,i)=Aa(:,4+mo*(i-1));
    gH(:,i)=Aa(:,5+mo*(i-1));
    gZ(:,i)=Aa(:,6+mo*(i-1));
    gCp(:,i)=Aa(:,7+mo*(i-1));
    gT(:,i)=Aa(:,8+mo*(i-1));
    gNTG(:,i)=Aa(:,9+mo*(i-1));
    gKY(:,i)=Aa(:,10+mo*(i-1));
    gKZ(:,i)=Aa(:,11+mo*(i-1));
end;
% gK=gK.*(gK>0)+0.1*(gK<0)+10;
% gK=gK/10;
v1=0:nl-1;
Won=repmat(1:size(WXY,1),nl,1)+repmat(v1'*size(XY,1),1,size(WXY,1)); 
Won=Won';
Won=Won(:);

BndZ=zeros(size(gKX(:),1),1);
BndZ(1:ncg)=1;
BndZ(ncg*(nl-1)+1:end)=2;

dH=gZ+gH;
NL=meshgrid(1:nl,1:ncg);
NamXY=repmat([1:ncg]',1,nl);

gKX(gKX<=0)=0.01;
gKY(gKY<=0)=0.01;
gKZ(gKZ<=0)=0.01;
gMp(gMp<=0)=0.001;
gSw(gSw<0)=0.01;
gH(gH<=0)=0.01;

DATA.gKX=gKX;
DATA.gKY=gKY;
DATA.gKZ=gKZ;
DATA.gMp=gMp;
DATA.gP=gP;
DATA.gSw=gSw;
DATA.gCp=gCp;
DATA.XY=XY;
DATA.gH=gH;
DATA.gZ=gZ;
DATA.gT=gT;
DATA.gNTG=gNTG;
DATA.Won=Won;
DATA.BndXY=BndXY;
DATA.BndZ=BndZ;
DATA.dH=dH;
DATA.NL=NL;
DATA.NamXY=NamXY;

end

function B=razmaz(A,WXY,XY)
 F=scatteredInterpolant(WXY(:,1),WXY(:,2),A);
 B=F(XY(:,1),XY(:,2));
%B=mean(A)*ones(size(XY,1),1);
end

function XY=NODuble(XY)
i=0;
while i<size(XY,1)
    i=i+1;
    a1=(XY(i,1)==XY(:,1)).*(XY(i,2)==XY(:,2));
    if sum(a1)>1
     r=find(a1);
     XY(r(2:end),:)=[];
    end;
end;
end