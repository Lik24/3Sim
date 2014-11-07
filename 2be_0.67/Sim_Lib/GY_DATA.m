function [GY_Data,KX,Sw,B,A2B,dVb,p]=GY_DATA(DATA,xy1,xy2,PR)

dh=24;
drob=6.5*dh;        % Густота сетки

bnd_XY=DATA.BndXY;
bnd_Z=DATA.BndZ;
P=DATA.gP;
KX=DATA.gKX;
Sw=DATA.gSw;

nz=size(bnd_Z);
nxy=size(bnd_XY);

GY_Kz=0*1.0*ones(nz)*8.64;
GY_Pz=39*ones(nz);
GY_Swz=1*ones(nz);

GY_Kz(bnd_Z~=2)=0;
GY_Kz(bnd_Z~=1)=0.628*ones(nxy)*8.64/20*0;

GY_Kxy=0.628*ones(nxy)*8.64*0.1*0;
GY_Pxy=P.*ones(nxy);
GY_Swxy=1*ones(nxy);

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
KX(IN==1)=mean(GY_Kxy);
Sw(IN==1)=mean(GY_Swxy);

[XY,GR] = Mesh4(xy1,drob,xy2);

% [IN1,ON1]=inpolygon(XY(:,1),XY(:,2),xy1(:,1),xy1(:,2));
[IN2,ON2]=inpolygon(XY(:,1),XY(:,2),xy2(:,1),xy2(:,2));
xy2=XY(ON2,:);
%%
[IN1,ON1]=inpolygon(XY(:,1),XY(:,2),xy1(:,1),xy1(:,2));
Con=XY(ON1==1,:);

[t,R]=cart2pol(Con(:,1)-500,Con(:,2)-1000);
[P,n]=sort(t);
XYgy1=Con(n,:);
        
for i=1:size(XYgy1,1)
 r(i,1)=find((XYgy1(i,1)==XY(:,1)).*(XYgy1(i,2)==XY(:,2)));
end;

BND1(:,1)=r;
BND1(:,2)=[r(2:end);r(1)];
r=[];
%%
[IN2,ON2]=inpolygon(XY(:,1),XY(:,2),xy2(:,1),xy2(:,2));
Con=XY(ON2==1,:);

[t,R]=cart2pol(Con(:,1)-500,Con(:,2)-1000);
[P,n]=sort(t);
XYgy2=Con(n,:);
        
for i=1:size(XYgy2,1)
 r(i,1)=find((XYgy2(i,1)==XY(:,1)).*(XYgy2(i,2)==XY(:,2)));
end;

BND2(:,1)=r;
BND2(:,2)=[r(2:end);r(1)];

BND=[BND1;BND2];

[B]=MR_Prop_Bond(XY,1,BND);
nb=size(B,1);

for i=1:size(XY,1)
  v1=XY(i,1)==DATA.XY(:,1);
  v2=XY(i,2)==DATA.XY(:,2);
  if sum(v1.*v2)~=0
      r1(i)=find(v1.*v2);
      c1(i)=i;
  end    
end;

Z=razmaz(DATA.gZ,DATA.XY,XY);
H=razmaz(DATA.gH,DATA.XY,XY);
K=razmaz(DATA.gKX,DATA.XY,XY);
P0=razmaz(DATA.gP,DATA.XY,XY);

[L,dB,S,H1]=Geome3_1(B,XY,Z,H);
r2=find(ON1);
B_GY=zeros(nb,1);

for i=r2'
  b=B(i,:);
  r3=find(b==1);
  on_b=ON1(r3);
  
  xy=XY(r3(on_b==1),:);
  B_GY(i,1)=sum(((xy(:,1)-XY(i,1)).^2+(xy(:,2)-XY(i,2)).^2).^0.5);
end


c2=find(ON2);
B_GY_2=zeros(nb,1);
for i=c2'
  b=B(i,:);
  r3=find(b==1);
  on_b=ON2(r3);
  
  xy=XY(r3(on_b==1),:);
  B_GY_2(i,1)=sum(((xy(:,1)-XY(i,1)).^2+(xy(:,2)-XY(i,2)).^2).^0.5);
end

H_B=H(c1).*B_GY(c1).*K(c1);
A2B=sparse(r1,c1,H_B,size(DATA.XY,1),size(B,1));
p=symrcm(B);
B=B(p,p);
L=L(p,p);
dB=dB(p,p);
S=S(p,p);
H1=H1(p,p);

dVb=full(sum(S,2).*H(p));

A2B=A2B(:,p);

[r,c]=find(B==1);  
KcKl=K(r,1)+K(c,1);
Ke=2*K(r,1).*K(c,1)./KcKl;

Lm=L(r+(c-1)*nb);
dBm=dB(r+(c-1)*nb);
H1m=H1(r+(c-1)*nb);
TB=Ke.*dBm.*H1m./Lm;
TB=TB.*PR.aw(1)./PR.mu(1);
TB=sparse(r,c,TB,nb,nb);
B=TB-sparse(1:nb,1:nb,sum(TB,2),nb,nb);

GY_Data.P0=P0;
GY_Data.T0=40;
GY_Data.XY=XY;
GY_Data.BND=BND1;

GY_Data.GY_Txy=H.*B_GY_2.*K;
end

function B=razmaz(A,WXY,XY)
 F=scatteredInterpolant(WXY(:,1),WXY(:,2),A,'nearest','nearest');
 B=F(XY(:,1),XY(:,2));
end
