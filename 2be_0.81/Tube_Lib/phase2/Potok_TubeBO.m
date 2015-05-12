function [CW,CO,CG,CORS,CP]=Potok_TubeBO(C,P,KWOG,Cp,PR,fp,kms,L,r,c,n,CMP,v)

Kwc=KWOG.w(v(c));
Kwl=KWOG.w(v(r));

Koc=KWOG.o(v(c));
Kol=KWOG.o(v(r));

Kgc=KWOG.g(v(c));
Kgl=KWOG.g(v(r));

Bw1 = CMP.Bw(v,2);
Bo1 = CMP.Bo(v,2);
Bg1 = CMP.Bg(v,2);
Rs1 = CMP.Rs(v,2);

vPw = P(c,1) - P(r,1)>0;
vPo = P(c,2) - P(r,2)>0;
vPg = P(c,3) - P(r,3)>0;

Cpc=Cp(c);
Cpl=Cp(r);

Cpe=Cpc.*vPw(r)+Cpl.*(vPw(r)==0);
Kfw=Kwc.*vPw(r)+Kwl.*(vPw(r)==0);
Kfo=Koc.*vPo(r)+Kol.*(vPo(r)==0);
Kfg=Kgc.*vPg(r)+Kgl.*(vPg(r)==0);
Rs = Rs1(c).*vPo(r) + Rs1(r).*(vPo(r)==0);

Bw = 0.5*(Bw1(c) + Bw1(r));
Bo = 0.5*(Bo1(c) + Bo1(r));
Bg = 0.5*(Bg1(c) + Bg1(r));

[CP,wmu]=poly_vis(Cpe,PR.mup,fp,Kfw,C,n,r,c,PR.mu(1));

% Tw=C.*Kfw.*((1-Cpe)./mu(1)+Cpe./mu(4));
Tw = C.*Kfw./Bw./wmu;
To = C.*Kfo./Bo/PR.mu(2);
Tg = C.*Kfg./Bg/PR.mu(3);      %//////////
%Tp=C.*Kfp./mu(4);
if kms~=0
     dPL=abs(dP./L((r+(c-1)*n)));
    %     проверить првильность задания плотности
    [Tw] = Forh(Tw,kms, PR.Ro(1), Kfw, PR.Kc, PR.mu(1), dPL);
    [To] = Forh(To,kms, PR.Ro(2), Kfo, PR.Kc, PR.mu(2), dPL);
end;

TW = sparse(r,c,Tw,n,n);       TW = TW + TW';
TORS = sparse(r,c,To.*Rs,n,n); TORS = TORS + TORS';
TO = sparse(r,c,To,n,n);       TO = TO + TO';
TG = sparse(r,c,Tg,n,n);       TG = TG + TG'; 

CW = TW - sparse(1:n,1:n,sum(TW,2),n,n);
CO = TO - sparse(1:n,1:n,sum(TO,2),n,n);
CORS = TORS - sparse(1:n,1:n,sum(TORS,2),n,n);
CG = TG - sparse(1:n,1:n,sum(TG,2),n,n);

