function [LA,LA0,LA_CND]=Potok_tepl(LT,Sw,La,ntg,mp,r,c,TL,TW,P,Cp,Ro,ndt,dt)
n=size(Sw,1);

% r=rc(:,1);
% c=rc(:,2);

Lfw=ntg.*mp.*Sw*La(1)+ntg.*mp.*(1-Sw)*La(2)+ntg.*(1-mp)*La(4)+(1-ntg)*La(5);
Lwc=Lfw(r);
Lwl=Lfw(c);

LT=LT.*2.*Lwc.*Lwl./(Lwc+Lwl);
LT(isnan(LT)==1)=0;

T1=sparse(r,c,LT,n,n);
LT=reshape(T1,n,n);

LA_CND=LT-sparse(1:n,1:n,sum(LT),n,n);

TO=TL-TW;
dP=P(r)-P(c);

Qw=TW(r+(c-1)*n).*dP*Cp(1)*Ro(1);
Qo=TO(r+(c-1)*n).*dP*Cp(2)*Ro(2);
% 
Qw1=TW(r+(c-1)*n).*dP.*Cp(1)*Ro(1).*(dP>=0);%(dP>=0);%*
Qo1=TO(r+(c-1)*n).*dP.*Cp(2)*Ro(2).*(dP>=0);%(dP>=0);%*

Qw2=TW(r+(c-1)*n).*dP.*Cp(1)*Ro(1).*(dP<0);%(dP<0);%*
Qo2=TO(r+(c-1)*n).*dP.*Cp(2)*Ro(2).*(dP<0);%(dP<0);%*

LT=sparse(r,c,Qo+Qw,n,n);
LT1=sparse(r,c,Qo1+Qw1,n,n);
LT2=sparse(r,c,Qo2+Qw2,n,n);

% sum(LT,2)
% sum(LT,1)
% kjgh
%sum(LT,2)
LA_cnv=LT1-sparse(1:n,1:n,sum(LT,2),n,n);
LA0=sum(LT2,1)*dt/ndt;

%LA0=0*LT1+sparse(1:n,1:n,sum(LT2,1),n,n);
% size(LA0)
% LA0(:,410)
% LT2(:,410)
LA=LA_cnv*dt/ndt;

%  full(LA)
% tyj