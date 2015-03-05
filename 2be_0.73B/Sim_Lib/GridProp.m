function [DATA]=GridProp(wKx,wKy,wKz,wMp,wP,wSw,wCp,wT,wNTG,CXY,wH,wZ,gXY,nl,WXY,XYgy,GY_subl,Won3)
[nw]=size(wKx,1);
A=zeros(nw,nl*6);
XY=[gXY];

XY=NODuble(XY);
%plot(XY(:,1),XY(:,2),'*')
% plot(gXY(:,1),gXY(:,2),'*')
% dfdf
[IN,ON]=inpolygon(XY(:,1),XY(:,2),XYgy(:,1),XYgy(:,2));
Con=XY(ON==1,:);

 xc=sum(Con(:,1),1)/size(Con,1);
 yc=sum(Con(:,2),1)/size(Con,1)+2000;
[t,R]=cart2pol(Con(:,1)-xc,Con(:,2)-yc);
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
DT = delaunayTriangulation(XY(:,1),XY(:,2),BND);
132 
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
[~,~,won1]=intersect(WXY,XY,'rows','stable');
Won(:,1)=repmat(won1,nl,1)+reshape(repmat(v1'*size(XY,1),1,size(WXY,1))',size(Won3,1)*nl,1); 
Won(:,2)=repmat(Won3(:,2),nl,1);
Won(:,3)=repmat(Won3(:,1),nl,1);
Won(:,4)=repmat(Won3(:,3),nl,1);

BndZ=zeros(size(gKX(:),1),1);
BndZ(1:ncg)=1;
BndZ(ncg*(nl-1)+1:end)=2;
%BndZ(XY(:,1)>1200)=0;

dH=gZ+gH;
NL=meshgrid(1:nl,1:ncg);
NamXY=repmat([1:ncg]',1,nl);

%GY=load('GY_Ura');
a=find(GY_subl(:,2)==-1);
for l=1:size(a,1)
 if l<size(a,1)   
 x_y(l)={GY_subl(a(l)+1:a(l+1)-1,:)};
 else
 x_y(l)={GY_subl(a(l)+1:end,:)};    
 end
end
% gKX(gKX<0)=0;
% gKY(gKY<0)=0;
% gKZ(gKZ<0)=0;
% gMp(gMp<0)=0;
% gSw(gSw<0)=0;
% gH(gH<=0.001)=0;
% gNTG(gNTG<=0)=0;

for l=1:size(gKX,2)
%  xy=GY.xy{l};
%  a=GY.a{l};

 gKX(:,l)=sub_bond(gKX(:,l),x_y{l},XY,1);
 gKY(:,l)=sub_bond(gKY(:,l),x_y{l},XY,1);
 gKZ(:,l)=sub_bond(gKZ(:,l),x_y{l},XY,1);
 gMp(:,l)=sub_bond(gMp(:,l),x_y{l},XY,1);
 gSw(:,l)=sub_bond(gSw(:,l),x_y{l},XY,0);
 gH(:,l)=sub_bond(gH(:,l),x_y{l},XY,1);
 gNTG(:,l)=sub_bond(gNTG(:,l),x_y{l},XY,1);
end

ka=(gKX+gKY+gKZ).*gMp.*gNTG.*gH~=0;

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
%dK=rand(size(gKX))-1;
%load('dk.mat','dK')

DATA.gKX=gKX;%+dK.*max(gKX)*0.5;
DATA.gKY=gKY;%+dK.*max(gKY)*0.5;
DATA.gKZ=gKZ;%+dK.*max(gKZ)*0.5;
DATA.gMp=gMp;
DATA.gP=gP+0.001*rand(size(gP));
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

function B=sub_bond(B,xy,XY,fl)

 in=inpolygon(XY(:,1),XY(:,2),xy(:,1),xy(:,2));
 tempKX=B(in==1);
 if fl==1
  tempKX(tempKX<=0)=mean(tempKX(tempKX>0));
 else
  tempKX(tempKX<0)=mean(tempKX(tempKX>=0));
 end
 B(in==1)=tempKX;
  
end