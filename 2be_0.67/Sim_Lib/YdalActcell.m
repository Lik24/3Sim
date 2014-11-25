function [A,L,B,S,H1,HV,XY,KX,KY,KZ,Mp,Sw,H,Z,P,MCp,Won,Wf]=YdalActcell(A,L,B,S,H1,HV,XY,KX,KY,KZ,Mp,Sw,H,Z,P,MCp,Won,Wf,ka)

A=A(ka==1,ka==1);
L=L(ka==1,ka==1);
B=B(ka==1,ka==1);
S=S(ka==1);
H1=H1(ka==1,ka==1);
HV=HV(ka==1);

XY=XY(ka==1);
KX=KX(ka==1);
KY=KY(ka==1);
KZ=KZ(ka==1);

Mp=Mp(ka==1);
Sw=Sw(ka==1);
H=H(ka==1);
Z=Z(ka==1);
P=P(ka==1);
MCp=MCp(ka==1);

ka1=ka(Won);
Wf=Wf(ka1==1);

im=zeros(size(ka));
im(Won)=1;
im=im(ka==1);
Won=find(im(:));
