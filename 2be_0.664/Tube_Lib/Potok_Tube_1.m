function [CW,CP,CG]=Potok_Tube_1(TC,P,vP1,vP2,Sw,Cp,PR,r,c)

ts=PR.ts;
tw=PR.tw;
mu=PR.mu;
bet=PR.bet;
Ro=PR.Ro;
dt=PR.dt;

Kfw=Sat_tube(Sw,1,1,ts,tw); %water

% Swc=Sw(r);
% Swl=Sw(c);

Kwc=Kfw(r);
Kwl=Kfw(c);

Cpc=Cp(r);
Cpl=Cp(c);

%Swe=Swc.*vP+Swl.*(vP==0);
Cpe=Cpc.*vP1+Cpl.*vP2;
Kf=Kwc.*vP1+Kwl.*vP2;

%Kf=Sat_tube(Swe,1,1,ts,tw); %water

CW=TC.*Kf.*((1-Cpe)./mu(1)+Cpe./mu(4));
CP=TC.*Kf.*Cpe./mu(4);

if bet~=0
    dP=(P(r)-P(c));
    CW=Iter_Fort(CW,dP,mu(1),Ro(1),bet,TC,Kf,dt);
    CP=Iter_Fort(CP,dP,mu(4),Ro(1),bet,TC,Kf,dt);
end;
%CW=sparse((r+(c-1)*n),ones(size(r)),Tw,n*n,1);
% Tw=reshape(Tw1,n,n);
% CW=Tw-sparse(1:n,1:n,sum(Tw),n,n);
CG=1;
end

function T=Iter_Fort(T,dP,mu,Ro,bet,C,Kf,dt)
t2=T;
t0=T;
fl=0;
j=0;
while (fl==0)*(j<15)==1
    j=j+1;
    W=abs(T.*dP)*dt;
    F=1./(1+bet*Ro*Kf.*C.*W/mu);
    T=t0.*F;
    fl=all(abs((T(t0~=0)-t2(t0~=0))./T(t0~=0))<=5e-2);
    %   cb(j)=sum(abs(T));
    %  cv(j,1)=sum(abs((T(t0~=0)-t2(t0~=0))./T(t0~=0)));
    t2=T;
end;
    %abs((T(t0~=0)-t2(t0~=0))./T(t0~=0))
  %  full(cv)
 %   full(cb)
end