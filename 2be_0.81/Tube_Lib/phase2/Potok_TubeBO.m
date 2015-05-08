function [CW,CO,CG,CP]=Potok_TubeBO(C,P,Sw,Cp,PR,fp,kms,L,r,c,n,CMP)

PW = P(:,1);
PO = P(:,2);
PG = P(:,3);

ts=PR.ts;
tw=PR.tw;
Kc=PR.Kc;

Kfw=Sat_tube(Sw,1,1,ts,tw); %water
Kfo=Sat_tube(Sw,2,1,ts,tw); %oil
Kfg=Sat_tube(Sw,3,1,ts,tw); %gas

vPW = PW(c) - PW(r)>0; 
vPO = PO(c) - PO(r)>0;
vPG = PG(c) - PG(r)>0;

Kwc=Kfw(c);
Kwl=Kfw(r);

Koc=Kfo(c);
Kol=Kfo(r);

Kgc=Kfg(c);
Kgl=Kfg(r);

Cpc=Cp(c);
Cpl=Cp(r);

Cpe=Cpc.*vPW+Cpl.*(vPW==0);
Kfw=Kwc.*vPW+Kwl.*(vPW==0);
Kfo=Koc.*vPO+Kol.*(vPO==0);
Kfg=Kgc.*vPO+Kgl.*(vPO==0);

Bw1 = CMP.Bw(:,2);
Bo1 = CMP.Bo(:,2);
Bg1 = CMP.Bg(:,2);
Rs1 = CMP.Rs(:,2);

Bw = 0.5*(Bw1(c) + Bw1(r));
Bo = 0.5*(Bo1(c) + Bo1(r));
Bg = 0.5*(Bg1(c) + Bg1(r));
Rs = 0.5*(Rs1(c) + Rs1(r));

[CP,wmu]=poly_vis(Cpe,PR.mup,fp,Kfw,C,n,r,c,PR.mu(1));

% Tw=C.*Kfw.*((1-Cpe)./mu(1)+Cpe./mu(4));
Tw = C.*Kfw./Bw./wmu;
To = C.*Kfo./Bo/PR.mu(2);
Tg = C.*Kfg./Bg/PR.mu(3);      %//////////
%Tp=C.*Kfp./mu(4);
if kms~=0
     dPL=abs(dP./L((r+(c-1)*n)));
    %     проверить првильность задания плотности
    [Tw] = Forh(Tw,kms, PR.Ro(1), Kfw, Kc, PR.mu(1), dPL);
    [To] = Forh(To,kms, PR.Ro(2), Kfo, Kc, PR.mu(2), dPL);
end;

Tw = sparse(r,c,Tw,n,n);       Tw = Tw + Tw';
To = sparse(r,c,To,n,n);       To = To + To';
Tg = sparse(r,c,Tg,n,n);       Tg = Tg + Tg'; 

CW = Tw - sparse(1:n,1:n,sum(Tw,2),n,n);
CO = To - sparse(1:n,1:n,sum(To,2),n,n);
CG = Tg - sparse(1:n,1:n,sum(Tg,2),n,n);

