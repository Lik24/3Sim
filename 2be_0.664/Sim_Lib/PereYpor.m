function [A,L,S,B,H1,K,XY,Mp,Sw,H,Z,P,MCp,T,NTG,p,rz,cz,BndXY,BndZ,dH,NL,NamXY,GYData1]=PereYpor(A,L,S,B,H1,KX,KY,KZ,Mp,Sw,XY,H,Z,P,MCp,DATA,GYData)

T=DATA.gT;
NTG=DATA.gNTG;
BndXY=DATA.BndXY;
BndZ=DATA.BndZ;
dH=DATA.dH;
NL=DATA.NL;
NamXY=DATA.NamXY;

A=A.*(A>0);

p=symrcm(A);
[r1,c1]=find(L);
[n,Nl]=size(KX);

A=A(p,p);
L=L(p,p);
S=S(p,p);
B=B(p,p); 
H1=H1(p,p);

KX=KX(:);
KY=KY(:);
KZ=KZ(:);
Mp=Mp(:);
Sw=Sw(:);
H=H(:); 
Z=Z(:);
P=P(:);
MCp=MCp(:);
T=T(:);
NTG=NTG(:);
dH=dH(:);
NamXY=NamXY(:);
NL=NL(:);

K(:,1)=KX(p);
K(:,2)=KY(p);
K(:,3)=KZ(p);

Mp=Mp(p);
Sw=Sw(p);
XY=XY(p,:);
H=H(p); 
Z=Z(p);
P=P(p); 
MCp=MCp(p);
T=T(p);
NTG=NTG(p);

dH=dH(p);
NamXY=NamXY(p);
NL=NL(p);

vr=[n+1:n*Nl,1:n*(Nl-1)];
vc=[1:n*(Nl-1),n+1:n*Nl];

AF=sparse(vr,vc,ones(1,2*n*(Nl-1)),n*Nl,n*Nl);
AZ=AF(p,p);
[rz,cz]=find(AZ);


GY_Kz=GYData.GY_Kz;
GY_Pz=GYData.GY_Pz;
GY_Swz=GYData.GY_Swz;

GY_Kxy=GYData.GY_Kxy;
GY_Pxy=GYData.GY_Pxy;
GY_Swxy=GYData.GY_Swxy;


BndXY=BndXY(p);
BndZ=BndZ(p);

GYData1.GY_Kz=GY_Kz(p);
GYData1.GY_Pz=GY_Pz(p);
GYData1.GY_Swz=GY_Swz(p);

GYData1.GY_Kxy=GY_Kxy(p);
GYData1.GY_Pxy=GY_Pxy(p);
GYData1.GY_Swxy=GY_Swxy(p);