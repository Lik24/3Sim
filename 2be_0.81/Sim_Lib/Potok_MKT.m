function [TW,TO,TG,TORS,TP]=Potok_MKT(T,P,KWOG,Cp,PR,rc,fp,kms,L,SGM,CMP)
n = size(P,1);

r = rc(:,1);
c = rc(:,2);

PW = P(:,1);
PO = P(:,2);
PG = P(:,3);

%vP = dP>0;
vPW = PW(c) - PW(r)>0; 
vPO = PO(c) - PO(r)>0;
vPG = PG(c) - PG(r)>0;

Kwc = KWOG.w(c);
Kwl = KWOG.w(r);

Koc = KWOG.o(c);
Kol = KWOG.o(r);

Kgc = KWOG.g(c);     %/////////////
Kgl = KWOG.g(r);  

Cpc = Cp(c);
Cpl = Cp(r);

%Swe=Swc.*vP+Swl.*(vP==0);
Cpe = Cpc.*vPW + Cpl.*(vPW==0);
Kfw1 = Kwc.*vPW + Kwl.*(vPW==0);
Kfo1 = Koc.*vPO + Kol.*(vPO==0);
Kfg1 = Kgc.*vPG + Kgl.*(vPG==0);        %///////

Bw = 0.5*(CMP.Bw(c,2) + CMP.Bw(r,2));
Bo = 0.5*(CMP.Bo(c,2) + CMP.Bo(r,2));
Bg = 0.5*(CMP.Bg(c,2) + CMP.Bg(r,2));
Rs = 0.5*(CMP.Rs(c,2) + CMP.Rs(r,2));

% Kfw=Sat_cal(Swe,1,1,as,aw); %water
% Kfo=Sat_cal(Swe,2,1,as,aw); %oil
% Kfp=Sat_cal(Swe.*Cpe,1,1,as,aw); %polim

[TP,wmu]=poly_vis(Cpe,PR.mup,fp,Kfw1,T,n,r,c,PR.mu(1));

% Tw=T.*Kfw1./mu(1);
Tw = T.*Kfw1./Bw./wmu;
To = T.*Kfo1./Bo/PR.mu(2);
Tg = T.*Kfg1./Bg/PR.mu(3);

 if kms~=0
     dPL=abs(dP./L((r+(c-1)*n)));
    %  проверить првильность задания плотности
    [Tw] = Forh(Tw,kms, PR.Ro(1), Kfw1, Ke, PR.mu(1), dPL);
    [To] = Forh(To,kms, PR.Ro(2), Kfo1, Ke, PR.mu(2), dPL);
end;

% Tw1=sparse((r+(c-1)*n),ones(size(r)),Tw,n*n,1);
% To1=sparse((r+(c-1)*n),ones(size(r)),To,n*n,1);


Tw = sparse(r,c,Tw,n,n);    Tw = Tw + Tw';
ToRs = sparse(r,c,To.*Rs,n,n);     ToRs = ToRs + ToRs';
%Trs = sparse(r,c,0.5*PR.rs.*To.*(PO(r)-PO(c)),n,n); Trs = Trs - Trs';
To = sparse(r,c,To,n,n);       To = To + To';
Tg = sparse(r,c,Tg,n,n);       Tg = Tg + Tg';     %//////////

TW = Tw - sparse(1:n,1:n,sum(Tw,2),n,n);
TO = To - sparse(1:n,1:n,sum(To,2),n,n);
TORS = ToRs - sparse(1:n,1:n,sum(ToRs,2),n,n);
%Trs = Trs + sparse(1:n,1:n,sum(Trs,2),n,n);
TG = Tg - sparse(1:n,1:n,sum(Tg,2),n,n);

% TP=Tp-sparse(1:n,1:n,sum(Tp,2),n,n);
% Tp1=sparse(r,c,Tp,n,n);
