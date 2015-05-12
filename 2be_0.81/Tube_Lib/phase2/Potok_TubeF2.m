function [CW,CO,CP]=Potok_TubeF2(TC,vP1,Kfw,Kfo,Cp,PR,r,c,kms,dP,L,n,CMP,v)

if isempty(TC)==0
mu=PR.mu;
Ro=PR.Ro;
Kc=PR.Kc;

Kwc=Kfw(r);
Kwl=Kfw(c);

Koc=Kfo(r);
Kol=Kfo(c);

Cpc=Cp(r);
Cpl=Cp(c);

%Swe=Swc.*vP+Swl.*(vP==0);
%sum([Cpc,Cpl].*[vP1,vP2]);

Cpe=Cpc.*vP1(:,1)+Cpl.*vP1(:,2);
Kfw=Kwc.*vP1(:,1)+Kwl.*vP1(:,2);
Kfo=Koc.*vP1(:,3)+Kol.*vP1(:,4);

Bw = 0.5*(CMP.Bw(v(c),2) + CMP.Bw(v(r),2));
Bo = 0.5*(CMP.Bo(v(c),2) + CMP.Bo(v(r),2));

%Kf=Sat_tube(Swe,1,1,ts,tw); %water

CO=TC.*Kfo./Bo./mu(2);
CW=TC.*Kfw./Bw.*((1-Cpe)./mu(1)+Cpe./mu(4));
CP=TC.*Kfw.*Cpe./mu(4);

if kms~=0
     dPL=abs(dP./L((r+(c-1)*n)));
    %     проверить првильность задания плотности
    [CW] = Forh(CW,kms, Ro(1), Kf, Kc, mu(1), dPL);
    [CO] = Forh(CO,kms, Ro(2), Kfo, Kc, mu(2), dPL);
end;


CW = sparse(r,c,CW,n,n);      CW = CW + CW';
CO = sparse(r,c,CO,n,n);       CO = CO + CO';

CW = CW - sparse(1:n,1:n,sum(CW,1),n,n);
CO = CO - sparse(1:n,1:n,sum(CO,1),n,n);
else
  CL=[];
  CW=[];
  CO=[];
  CP=[];
end;
end

function T=Iter_Fort(T,dP,mu,Ro,bet,C,Kf,dt)
t2=T;
t0=T;
fl=1>0;
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