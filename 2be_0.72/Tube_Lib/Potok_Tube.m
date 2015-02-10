function [CL,CW,CP,CG]=Potok_Tube(C,P,Sw,Cp,PR,mup,fp,kms,L,r,c,n,A,DZ)

ts=PR.ts;
tw=PR.tw;
mu=PR.mu;
Ro=PR.Ro;
Kc=PR.Kc;
dZW=DZ{1};  dZO=DZ{2};  dZG=DZ{3};

Kfw=Sat_tube(Sw,1,1,ts,tw); %water
Kfo=Sat_tube(Sw,2,1,ts,tw); %oil

dP=P(c)-P(r);
vP=dP>0;

vPW=dP+dZW(r+(c-1)*n)>0;
vPO=dP+dZO(r+(c-1)*n)>0;
vPG=dP+dZG(r+(c-1)*n)>0;

Kwc=Kfw(c);
Kwl=Kfw(r);

Koc=Kfo(c);
Kol=Kfo(r);

Cpc=Cp(c);
Cpl=Cp(r);

Cpe=Cpc.*vPW+Cpl.*(vPW==0);
Kfw=Kwc.*vPW+Kwl.*(vPW==0);
Kfo=Koc.*vPO+Kol.*(vPO==0);

[CP,wmu]=poly_vis(Cpe,mup,fp,Kfw,C,n,r,c,mu(1));


% Tw=C.*Kfw.*((1-Cpe)./mu(1)+Cpe./mu(4));
Tw=C.*Kfw./wmu;
To=C.*Kfo./mu(2);

%Tp=C.*Kfp./mu(4);
if kms~=0
     dPL=abs(dP./L((r+(c-1)*n)));
    %     проверить првильность задания плотности
    [Tw] = Forh(Tw,kms, Ro(1), Kfw, Kc, mu(1), dPL);
    [To] = Forh(To,kms, Ro(2), Kfo, Kc, mu(2), dPL);
end;

Twa=Tw.*A(c);
Tw=sparse(r,c,Tw,n,n);  Tw=Tw+Tw';
Twa=sparse(r,c,Twa,n,n);  Twa=Twa+Twa';
To=sparse(r,c,To,n,n);  To=To+To';

CL=To+Twa-sparse(1:n,1:n,sum(To+Twa,2),n,n);
CW=Tw-sparse(1:n,1:n,sum(Tw,2),n,n);
CG=1;
