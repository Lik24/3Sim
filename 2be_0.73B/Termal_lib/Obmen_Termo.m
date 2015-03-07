function [A2CLi,A2CLo]=Obmen_Termo(A2C,Sw,Cw,La,ntg,mp,mc,r,c,na,nc,Pm,Pt,TL,TW,Cp,Ro,ndt,dt)

Lfa=ntg.*mp.*Sw*La(1)+ntg.*mp.*(1-Sw)*La(2)+ntg.*(1-mp)*La(4)+(1-ntg)*La(5);
Lfc=mc.*Cw*La(1)+mc.*(1-Cw)*La(2)+(1-mc)*La(4);

L=A2C(r).*(Lfa(r)+Lfc(c))/2;
A2CL=sparse(r,c,L,na,nc);

[r,c]=find(TL);
dP=Pm(r)-Pt(c);

TO=TL-TW;
Qw=TW(r+(c-1)*na).*dP*Cp(1)*Ro(1).*(dP>0);
Qo=TO(r+(c-1)*na).*dP*Cp(2)*Ro(2).*(dP>0);

Qw1=TW(r+(c-1)*na).*dP*Cp(1)*Ro(1).*(dP<0);
Qo1=TO(r+(c-1)*na).*dP*Cp(2)*Ro(2).*(dP<0);

A2CLo=-sparse(r,c,Qw+Qo,na,nc)*dt/ndt;
A2CLi=-sparse(r,c,Qw1+Qo1,na,nc)*dt/ndt;
% TL
% [r,c,r+(c-1)*na]
%  full([sum(A2CLi,2),sum(A2CLo,2)])
%   fgh

A2CLi=A2CLi+A2CL;
A2CLo=A2CLo+A2CL;

