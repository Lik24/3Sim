function [CW,CO,CP]=Potok_Tube2(C,P,KWOG,Cp,PR,fp,kms,L,r,c,n,CMP,v)

PW = P(:,1);
PO = P(:,2);

ts=PR.ts;
tw=PR.tw;
Kc=PR.Kc;

vPW = PW(c) - PW(r)>0; 
vPO = PO(c) - PO(r)>0;

Kwc=KWOG.w(v(c));
Kwl=KWOG.w(v(r));

Koc=KWOG.o(v(c));
Kol=KWOG.o(v(r));

Cpc=Cp(c);
Cpl=Cp(r);

Cpe=Cpc.*vPW+Cpl.*(vPW==0);
Kfw=Kwc.*vPW+Kwl.*(vPW==0);
Kfo=Koc.*vPO+Kol.*(vPO==0);

Bw1 = CMP.Bw(:,2);
Bo1 = CMP.Bo(:,2);

Bw = 0.5*(Bw1(v(r)) + Bw1(v(c)));
Bo = 0.5*(Bo1(v(r)) + Bo1(v(c)));

[CP,wmu]=poly_vis(Cpe,PR.mup,fp,Kfw,C,n,r,c,PR.mu(1));

% Tw=C.*Kfw.*((1-Cpe)./mu(1)+Cpe./mu(4));
Tw = C.*Kfw./Bw./wmu;
To = C.*Kfo./Bo/PR.mu(2);
%Tp=C.*Kfp./mu(4);
if kms~=0
     dPL=abs(dP./L((r+(c-1)*n)));
    %     проверить првильность задания плотности
    [Tw] = Forh(Tw,kms, PR.Ro(1), Kfw, Kc, PR.mu(1), dPL);
    [To] = Forh(To,kms, PR.Ro(2), Kfo, Kc, PR.mu(2), dPL);
end;


Tw = sparse(r,c,Tw,n,n);       Tw = Tw + Tw';
To = sparse(r,c,To,n,n);       To = To + To';

CW = Tw - sparse(1:n,1:n,sum(Tw,2),n,n);
CO = To - sparse(1:n,1:n,sum(To,2),n,n);

