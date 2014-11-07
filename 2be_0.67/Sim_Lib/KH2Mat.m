function [Ke,Ke_gy_XY,Ke_gy_Z,dV]=KH2Mat(K,h,Mp,S,r,c,rz,cz,Kz_gy,Kxy_gy,bz,BndXY,BndZ)
% H=(h(r)+h(c))/2;
K(K==0)=0.01;
n=size(K,1);
KcKl=K(r,1)+K(c,1);
Ke=2*K(r,1).*K(c,1)./KcKl;

KcKl=K(rz,3)+K(cz,3);
Kez=2*K(rz,3).*K(cz,3)./KcKl;
% Ke(rz+(cz-1)*n)=Kez;
% size(Ke)
% size(r)
Ke(Ke==0)=-1;
KE=sparse(r,c,Ke,n,n);
%Kez'
Kez(Kez==0)=-1;
KE(rz+(cz-1)*n)=Kez;


[r,c,Ke]=find(KE);
Ke(isnan(Ke)==1)=0;
Ke(Ke==-1)=0;

dV=full(Mp.*sum(S,2).*h);

Ke_gy_Z=2*K(:,3).*Kz_gy./(K(:,3)+Kz_gy);
Ke_gy_Z(isnan(Ke_gy_Z)==1)=0;

Ke_gy_XY=2*K(:,1).*Kxy_gy./(K(:,1)+Kxy_gy);
Ke_gy_XY(isnan(Ke_gy_XY)==1)=0;
Ke_gy_XY=Ke_gy_XY.*BndXY;
