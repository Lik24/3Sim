function [TL,TW,TG]=Potok_MKT(T,P,Sw,Cp,as,aw,mu,rc)
n=size(P,1);

r=rc(:,1);
c=rc(:,2);

Kfw=Sat_cal(Sw,1,1,as,aw); %water
Kfo=Sat_cal(Sw,2,1,as,aw); %oil

vP=P(c)>P(r);
% Swc=Sw(r);
% Swl=Sw(c);

Kwc=Kfw(c);
Kwl=Kfw(r);

Koc=Kfo(c);
Kol=Kfo(r);

Cpc=Cp(c);
Cpl=Cp(r);

% sw=c.*vP+r.*(vP==0);
% SD=sparse(r,c,sw);
% P(4:5)
% [r,c,vP]
%SD
%sum(SD,1)-sum(SD,2)'

%Swe=Swc.*vP+Swl.*(vP==0);
%Cpe=Cpc.*vP+Cpl.*(vP==0);
Kfw1=Kwc.*vP+Kwl.*(vP==0);
Kfo1=Koc.*vP+Kol.*(vP==0);

% Kfw=Sat_cal(Swe,1,1,as,aw); %water
% Kfo=Sat_cal(Swe,2,1,as,aw); %oil
% Kfp=Sat_cal(Swe.*Cpe,1,1,as,aw); %polim

Tw=T.*Kfw1./mu(1);
To=T.*Kfo1./mu(2);
%Tp=T.*Kfp./mu(4);

% Tw1=sparse((r+(c-1)*n),ones(size(r)),Tw,n*n,1);
% To1=sparse((r+(c-1)*n),ones(size(r)),To,n*n,1);

Tw1=sparse(r,c,Tw,n,n);
To1=sparse(r,c,To,n,n);
%Tp1=sparse(r,c,Tp,n,n);

Tw=reshape(Tw1,n,n);
To=reshape(To1,n,n);
%Tp=reshape(Tp1,n,n);
% sum(To1+Tw1,2)-sum(To1+Tw1)'
% jkkljlk

TL=To+Tw-sparse(1:n,1:n,sum(To+Tw,2),n,n);
TW=Tw-sparse(1:n,1:n,sum(Tw,2),n,n);
TG=1;