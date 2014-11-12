function [DATA]=GridProp(wKx,wKy,wKz,wMp,wP,wSw,wCp,wT,wNTG,CXY,wH,wZ,gXY,nl,WXY,XYgy)
[nw]=size(wKx,1);
A=zeros(nw,nl*6);
XY=[gXY];

XY=NODuble(XY);
%plot(XY(:,1),XY(:,2),'*')
% plot(gXY(:,1),gXY(:,2),'*')
% dfdf
[IN,ON]=inpolygon(XY(:,1),XY(:,2),XYgy(:,1),XYgy(:,2));
Con=XY(ON==1,:);

% xc=sum(Con(:,1),1)/size(Con,1);
% yc=sum(Con(:,2),1)/size(Con,1);
[t,R]=cart2pol(Con(:,1)-15000,Con(:,2)-5000);
[P,n]=sort(t);
XYgy=Con(n,:);
        
% [IN,ON]=inpolygon(XY(:,1),XY(:,2),XYc(:,1),XYc(:,2));
% Con=XY(ON==1,:);
%         r=r(n) 
for i=1:size(XYgy,1)
 r(i,1)=find((XYgy(i,1)==XY(:,1)).*(XYgy(i,2)==XY(:,2)));
end;

BND(:,1)=r;
BND(:,2)=[r(2:end);r(1)];
%DT = delaunayTriangulation(XY(:,1),XY(:,2),BND);
 
% plot(XYgy(:,1),XYgy(:,2))
% plot(XY(BND,1),XY(BND,2))
% triplot(DT)

% XY=DT.Points;
% [IN,ON]=inpolygon(XY(:,1),XY(:,2),XYc(:,1),XYc(:,2));
% r=find(ON);
% BND1(:,1)=r;
% BND1(:,2)=[r(2:end);r(1)];
% 
% DT1 = delaunayTriangulation(XY(:,1),XY(:,2),BND1); 
%  IO = isInterior(DT1); 
% BND=BND1; 
%   triplot(DT1)
%    hold on
% for i=1:size(BND,1)
%   plot(XY(BND(:,:),1),XY(BND(:,:),2))
%    hold on
%    pause(0.1)
% end
%   kjjh
ncg=size(XY,1);

BndXY=zeros(ncg,1);
ao=r;
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
%BndZ(XY(:,1)>1200)=0;

dH=gZ+gH;
NL=meshgrid(1:nl,1:ncg);
NamXY=repmat([1:ncg]',1,nl);

gKX(gKX<=0)=0;
gKY(gKY<=0)=0;
gKZ(gKZ<=0)=0;
gMp(gMp<=0)=0;
gSw(gSw<0)=0;
gH(gH<=0.001)=0;
gNTG(gNTG<=0)=0;

ka=(gKX+gKY+gKZ).*gMp.*gNTG.*gH;

% gSw(XY(:,1)<500)=0.3;
% gSw(XY(:,2)>1150)=0.3;


% gKX(XY(:,1)<500)=mean(wKx);
% gKX(XY(:,2)>1150)=mean(wKx);
% 
% gMp(XY(:,1)<500)=0.3289;
% gMp(XY(:,2)>1150)=0.3289;

% gH(XY(:,1)<500)=10.94;
% gH(XY(:,2)>1150)=10.94;
% 
% gNTG(XY(:,1)<500)=6.62/10.94;
% gNTG(XY(:,2)>1150)=6.62/10.94;

% mean([gKX(XY(:,1)<500);gKX(XY(:,2)>1050)])
% mean([gMp(XY(:,1)<500);gMp(XY(:,2)>1050)])
% mean([gH(XY(:,1)<500);gH(XY(:,2)>1050)])
% mean([gNTG(XY(:,1)<500);gNTG(XY(:,2)>1050)])
%[gSw0,gSw]

DATA.gKX=gKX;
DATA.gKY=gKY;
DATA.gKZ=gKZ;
DATA.gMp=gMp;
DATA.gP=gP;
DATA.gSw=gSw;
DATA.gCp=gCp;
DATA.XY=XY;
DATA.BND=BND;
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
DATA.ka=ka>0;
DATA.XYgy=XYgy;

end

function B=razmaz(A,WXY,XY)
 F=scatteredInterpolant(WXY(:,1),WXY(:,2),A,'linear','nearest');
 B=F(XY(:,1),XY(:,2));
  
% mx(1)=min(WXY(:,1));
% mx(2)=max(WXY(:,1));
% my(1)=min(WXY(:,2));
% my(2)=max(WXY(:,2));
% [X,Y]=meshgrid(mx(1):50:mx(2),my(1):50:my(2));%,mz(1):5:mz(2)
% A(A==0)=nan;
% F=scatteredInterpolant(WXY(:,1),WXY(:,2),A,'linear','nearest');
% Vz=F(X,Y);
% Vz=inpaintn(Vz);
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