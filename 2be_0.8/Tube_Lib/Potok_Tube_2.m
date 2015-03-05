function [CL,CW,CO,CG,CWCgs,COCgs,CP]=Potok_Tube_2(TC,P,vP1,Kfw,Kfo,Kfg,Cp,PR,r,c,kms,dP,L,n,SGM)

if isempty(TC)==0
mu=PR.mu;
Ro=PR.Ro;
Kc=PR.Kc;

Kwc=Kfw(r);
Kwl=Kfw(c);

Koc=Kfo(r);
Kol=Kfo(c);

Kgc=Kfg(r);
Kgl=Kfg(c);

Cpc=Cp(r);
Cpl=Cp(c);

%Swe=Swc.*vP+Swl.*(vP==0);
%sum([Cpc,Cpl].*[vP1,vP2]);

Cpe=Cpc.*vP1(:,1)+Cpl.*vP1(:,2);
Kfw=Kwc.*vP1(:,1)+Kwl.*vP1(:,2);
Kfo=Koc.*vP1(:,3)+Kol.*vP1(:,4);
Kfg=Kgc.*vP1(:,5)+Kgl.*vP1(:,6);

Bw = 0.5*(SGM.Bwog(c,2) + SGM.Bwog(r,2));
Bo = 0.5*(SGM.Bwog(c,4) + SGM.Bwog(r,4));
Bg = 0.5*(SGM.Bwog(c,6) + SGM.Bwog(r,6));
Rs = 0.5*(SGM.Rs(c) + SGM.Rs(r));

%Kf=Sat_tube(Swe,1,1,ts,tw); %water

CO=TC.*Kfo./Bo./mu(2);
CG=TC.*Kfg./Bg./mu(3);
CW=TC.*Kfw./Bw.*((1-Cpe)./mu(1)+Cpe./mu(4));
CP=TC.*Kfw.*Cpe./mu(4);

if kms~=0
     dPL=abs(dP./L((r+(c-1)*n)));
    %     проверить првильность задания плотности
    [CW] = Forh(CW,kms, Ro(1), Kf, Kc, mu(1), dPL);
    [CO] = Forh(CO,kms, Ro(2), Kfo, Kc, mu(2), dPL);
end;

CWa = -CW.*SGM.Cgsw(c);
CWa = sparse(r,c,CWa,n,n);     CWa = CWa + CWa';
CW = sparse(r,c,CW,n,n);      CW = CW + CW';

COa = CO.*(Rs - SGM.Cgso(c));
COa = sparse(r,c,COa,n,n);     COa = COa + COa';
CO = sparse(r,c,CO,n,n);       CO = CO + CO';

CG = sparse(r,c,CG,n,n);       CG = CG + CG';     %//////////

CW = CW - sparse(1:n,1:n,sum(CW,1),n,n);
CO = CO - sparse(1:n,1:n,sum(CO,1),n,n);
CWCgs = CWa - sparse(1:n,1:n,sum(CWa,1),n,n);
COCgs = COa - sparse(1:n,1:n,sum(COa,1),n,n);
CG = CG - sparse(1:n,1:n,sum(CG,1),n,n);
CL = CG + CWCgs + COCgs;

else
  CL=[];
  CW=[];
  CO=[];
  CP=[];
  CG=[];
  COCgs=[];
  CWCgs=[];
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