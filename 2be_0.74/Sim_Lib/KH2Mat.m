function [Ke,Ke_gy,dV]=KH2Mat(K,H,Mp,S,r,c,rz,cz,Kz_gy,Kxy_gy,BndXY,DData,ka)
% H=(h(r)+h(c))/2;
%K(K==0)=0.01;
n=numel(K);
KcKl=K(r,1)+K(c,1);
Ke=2*K(r,1).*K(c,1)./KcKl;

KcKl=K(rz,3)+K(cz,3);
Kez=2*K(rz,3).*K(cz,3)./KcKl;

Ke(Ke==0)=-1;
KE=sparse(r,c,Ke,n,n);

%Kez(Kez==0)=-1;
KE(rz+(cz-1)*n)=Kez;

Ke=KE(r+(c-1)*n);
% 
% [r,c,Ke]=find(KE);
% Ke(isnan(Ke)==1)=0;
% Ke(Ke==-1)=0;

dV=full(Mp.*S.*H);

Ke_gy(:,2)=2*K(:,3).*Kz_gy./(K(:,3)+Kz_gy);
Ke_gy(isnan(Ke_gy(:,2))==1,2)=0;

Ke_gy(:,1)=2*K(:,1).*Kxy_gy./(K(:,1)+Kxy_gy);
Ke_gy(isnan(Ke_gy(:,1))==1,1)=0;
Ke_gy(:,1)=Ke_gy(:,1).*BndXY;


% Kd=DData.Kd;
% Kz_gy=Kz_gy.*Kd(:,1)./K(:,1);
% Kxy_gy=Kxy_gy.*Kd(:,3)./K(:,3);
% 
% Ke_gy(:,4)=2*Kd(:,3).*Kz_gy./(Kd(:,3)+Kz_gy);
% Ke_gy(isnan(Ke_gy(:,2))==1,4)=0;
% 
% Ke_gy(:,3)=2*Kd(:,1).*Kxy_gy./(Kd(:,1)+Kxy_gy);
% Ke_gy(isnan(Ke_gy(:,1))==1,3)=0;
% Ke_gy(:,1)=Ke_gy(:,1).*BndXY;