function [GY_Data,KX,Sw,B,A2B,dVb,p]=GY_DATA(on,DATA,xy2,PR)

drob=PR.drob; % Густота сетки

bnd_XY=DATA.BndXY;
bnd_Z=DATA.BndZ;
P=DATA.gP;
KX=DATA.gKX;
Sw=DATA.gSw;

nz=size(bnd_Z);
nxy=size(bnd_XY);
nxy(1)=nxy(1)/PR.Nl;

GY_Kz=1.0*ones(nz)*8.64;
GY_Pz=39*ones(nz);
GY_Swz= ones(nz);

GY_Kz(bnd_Z~=2)=0;
GY_Kz(bnd_Z~=1)=0.628*8.64/20*0;

GY_Kxy=ones(size(P))*8.64*0;
%GY_Kxy(repmat(DATA.XY(:,2),PR.Nl,1)>0)=0;
% GY_Kxy(DATA.XY(:,2)<1)=8.64*1;
% GY_Kxy(DATA.XY(:,1)>999)=8.64*1;

GY_Pxy=100*ones(size(P));%P.*ones(size(P));
GY_Swxy =  zeros(size(P));

GY_Data.GY_Kz=GY_Kz;
GY_Data.GY_Pz=GY_Pz;
GY_Data.GY_Swz=GY_Swz;

GY_Data.GY_Kxy=GY_Kxy;
GY_Data.GY_Pxy=GY_Pxy;
GY_Data.GY_Swxy=GY_Swxy;

xy1=DATA.XYgy;
xy1(end+1,:)=xy1(1,:);
[IN1]=inpolygon(DATA.XY(:,1),DATA.XY(:,2),xy1(:,1),xy1(:,2));
[IN2]=inpolygon(DATA.XY(:,1),DATA.XY(:,2),xy2(:,1),xy2(:,2));

IN=IN2-IN1;
for l=1:PR.Nl
KX(IN==1,l)=mean(GY_Kxy(:,l));
Sw(IN==1,l)=mean(GY_Swxy(:,l));
end

[XY,GR] = Mesh4(xy1,drob,xy2);

% [IN1,ON1]=inpolygon(XY(:,1),XY(:,2),xy1(:,1),xy1(:,2));
[IN2,ON2]=inpolygon(XY(:,1),XY(:,2),xy2(:,1),xy2(:,2));
xy2=XY(ON2,:);
%%
[~,ON1]=inpolygon(XY(:,1),XY(:,2),xy1(:,1),xy1(:,2));
Con=XY(ON1==1,:);

[t,R]=cart2pol(Con(:,1)-250,Con(:,2)-250);
[P,n]=sort(t);
XYgy1=Con(n,:);

r=[];
for i=1:size(XYgy1,1)
 r(i,1)=find((XYgy1(i,1)==XY(:,1)).*(XYgy1(i,2)==XY(:,2)));
end;
if isempty(r)==0
BND1(:,1)=r;
BND1(:,2)=[r(2:end);r(1)];
else
 BND1=XY;    
end
r=[];
%%
[IN2,ON2]=inpolygon(XY(:,1),XY(:,2),xy2(:,1),xy2(:,2));
Con=XY(ON2==1,:);

[t,R]=cart2pol(Con(:,1)-250,Con(:,2)-250);
[P,n]=sort(t);
XYgy2=Con(n,:);
        
for i=1:size(XYgy2,1)
 r(i,1)=find((XYgy2(i,1)==XY(:,1)).*(XYgy2(i,2)==XY(:,2)));
end;

if isempty(r)==0
BND2(:,1)=r;
BND2(:,2)=[r(2:end);r(1)];
else
 BND2=XY;    
end

BND=[BND1;BND2];

[B]=MR_Prop_Bond(XY,PR.Nl,BND);
nb=size(B,1);

% for i=1:size(XY,1)
%   v1=XY(i,1)==DATA.XY(:,1);
%   v2=XY(i,2)==DATA.XY(:,2);
%   if sum(v1.*v2)~=0
%       r1(i)=find(v1.*v2);
%       c1(i)=i;
%   end    
% end;
[~,c1,r1]=intersect(XY,DATA.XY,'rows','stable');

for l=1:PR.Nl
    Z(:,l)=razmaz(DATA.gZ(:,l),DATA.XY,XY);
    H(:,l)=razmaz(DATA.gH(:,l),DATA.XY,XY);
    K(:,l)=razmaz(DATA.gKX(:,l),DATA.XY,XY);
    P0(:,l)=razmaz(DATA.gP(:,l),DATA.XY,XY);
end

[L,dB,S,H1]=Geome3_1(B,XY,Z,H);
r2=find(ON1);

n=size(XY,1);
B1=B(1:n,1:n);
B_GY=zeros(nb/n,1);

for i=r2'
  b=B1(i,:);
  r3=find(b==1);
  on_b=ON1(r3);
  
  xy=XY(r3(on_b==1),:);
  B_GY(i,1)=sum(((xy(:,1)-XY(i,1)).^2+(xy(:,2)-XY(i,2)).^2).^0.5);
end


c2=find(ON2);
B_GY_2=zeros(nb/n,1);

for i=c2'
  b=B1(i,:);
  r3=find(b==1);
  on_b=ON2(r3);
  xy=XY(r3(on_b==1),:);
  B_GY_2(i,1)=sum(((xy(:,1)-XY(i,1)).^2+(xy(:,2)-XY(i,2)).^2).^0.5);
end

a2b={sparse(size(DATA.XY,1),size(B1,1))};
a2b=repmat(a2b,PR.Nl,PR.Nl);
for l=1:PR.Nl
 H_B=H(c1,l).*B_GY(c1).*K(c1,l);
 a2b(l,l)={sparse(r1,c1,H_B,size(DATA.XY,1),size(B1,1))};
end
A2B=cell2mat(a2b);

p=symrcm(B);
B=B(p,p);
L=L(p,p);
dB=dB(p,p);
S=S(p);
H1=H1(p,p);

Hp(:,1)=H(p);
dVb=full(S.*Hp);

A2B=A2B(:,p);

[r,c]=find(B==1);  
KcKl=K(r)+K(c);
Ke=2*K(r).*K(c)./KcKl;

Lm=L(r+(c-1)*nb);
dBm=dB(r+(c-1)*nb);
H1m=H1(r+(c-1)*nb);
TB=Ke.*dBm.*H1m./Lm;
TB=TB.*PR.aw(1)./PR.mu(1);
TB=sparse(r,c,TB,nb,nb);
B=TB-sparse(1:nb,1:nb,sum(TB,2),nb,nb);

GY_Data.P0=zeros(0,1);
GY_Data.T0=[];
GY_Data.XY=XY;
GY_Data.BND=[[],[]];

GY_Data.GY_Txy=zeros(0,1);%H.*K;%.*B_GY_2;%zeros(0,1);

if on==0
    B=[];
    A2B(:,1:end)=[];
    dVb=[];
    p=[];
end
end

function B=razmaz(A,WXY,XY)
 F=scatteredInterpolant(WXY(:,1),WXY(:,2),A,'nearest','nearest');
 B=F(XY(:,1),XY(:,2));
end
