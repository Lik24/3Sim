function [CL,CW,CG]=Potok_Tube(Cm,P,Sw,Cp,PR)

ts=PR.ts;
tw=PR.tw;
mu=PR.mu;
bet=PR.bet;
Ro=PR.Ro;

n=size(Cm,1);
[r,c]=find(Cm);
C=Cm(r+(c-1)*n);

Kfw=Sat_tube(Sw,1,1,ts,tw); %water
Kfo=Sat_tube(Sw,2,1,ts,tw); %oil

vP=P(r)>=P(c);
% Swc=Sw(r);
% Swl=Sw(c);

Kwc=Kfw(r);
Kwl=Kfw(c);

Koc=Kfo(r);
Kol=Kfo(c);

Cpc=Cp(r);
Cpl=Cp(c);

%Swe=Swc.*vP+Swl.*(vP==0);
Cpe=Cpc.*vP+Cpl.*(vP==0);
Kfw=Kwc.*vP+Kwl.*(vP==0);
Kfo=Koc.*vP+Kol.*(vP==0);

% Kfw=Sat_tube(Swe,1,1,ts,tw); %water
% Kfo=Sat_tube(Swe,2,1,ts,tw); %oil
%Kfp=Sat_tube(Swe.*Cpe,1,1,ts,tw); %polim

Tw=C.*Kfw.*((1-Cpe)./mu(1)+Cpe./mu(4));
To=C.*Kfo./mu(2);
%Tp=C.*Kfp./mu(4);

if bet~=0
    dP=(P(r)-P(c));
    Tw=Iter_Fort(Tw,dP,mu(1),Ro(1),bet,C,Kfw);
    To=Iter_Fort(To,dP,mu(2),Ro(2),bet,C,Kfo);
  %  Tp=Iter_Fort(Tp,dP,mu(4),Ro(1),bet,C,Kfp);
end;

Tw1=sparse(r,c,Tw,n,n);
To1=sparse(r,c,To,n,n);
%Tp1=sparse(r,c,Tp,n,n);

Tw=reshape(Tw1,n,n);
To=reshape(To1,n,n);
%Tp=reshape(Tp1,n,n);

CL=To+Tw-sparse(1:n,1:n,sum(To+Tw,2),n,n);
CW=Tw-sparse(1:n,1:n,sum(Tw,2),n,n);
CG=1;
end

function T=Iter_Fort(T,dP,mu,Ro,bet,C,Kf)
t2=T;
t0=T;
fl=0;
j=0;
while (fl==0)*(j<15)==1
    j=j+1;
    W=abs(T.*dP);
    F=1./(1+bet*Ro*Kf.*C.*W/mu);
    T=t0.*F;
    fl=all(abs((T(t0~=0)-t2(t0~=0))./T(t0~=0))<=5e-2);
  %   cb(j)=sum(abs(T));
  %  cv(j,1)=sum(abs((T(t0~=0)-t2(t0~=0))./T(t0~=0)));   
    t2=T;
end;

end