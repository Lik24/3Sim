function [CL,CW,CO,CG,CWCgs,COCgs,CP]=Potok_Tube(C,P,Sw,Cp,PR,fp,kms,L,r,c,n,SGM)

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

Bw = SGM.Bwog(:,2);
Bo = SGM.Bwog(:,4);
Bg = SGM.Bwog(:,6);
Rs = SGM.Rs(:,2);

Bw1 = 0.5*(Bw(c) + Bw(r));
Bo1 = 0.5*(Bo(c) + Bo(r));
Bg1 = 0.5*(Bg(c) + Bg(r));
Rs1 = 0.5*(Rs(c) + Rs(r));

[CP,wmu]=poly_vis(Cpe,PR.mup,fp,Kfw,C,n,r,c,PR.mu(1));

% Tw=C.*Kfw.*((1-Cpe)./mu(1)+Cpe./mu(4));
Tw = C.*Kfw./Bw1./wmu;
To = C.*Kfo./Bo1/PR.mu(2);
Tg = C.*Kfg./Bg1/PR.mu(3);      %//////////

%Tp=C.*Kfp./mu(4);
if kms~=0
     dPL=abs(dP./L((r+(c-1)*n)));
    %     проверить првильность задания плотности
    [Tw] = Forh(Tw,kms, PR.Ro(1), Kfw, Kc, PR.mu(1), dPL);
    [To] = Forh(To,kms, PR.Ro(2), Kfo, Kc, PR.mu(2), dPL);
end;

Twa = -Tw.*SGM.Cgsw(c);
Tw1a = sparse(r,c,Twa,n,n);     Tw1a = Tw1a + Tw1a';
Tw1 = sparse(r,c,Tw,n,n);       Tw1 = Tw1 + Tw1';

Toa = -To.*SGM.Cgso(c);
To1a = sparse(r,c,Toa,n,n);     To1a = To1a + To1a';
ToRs = To.*Rs1;
To1Rs = sparse(r,c,ToRs,n,n);   To1Rs = To1Rs + To1Rs';
To1 = sparse(r,c,To,n,n);       To1 = To1 + To1';

Tg1 = sparse(r,c,Tg,n,n);       Tg1 = Tg1 + Tg1';     %//////////

CW = Tw1 - sparse(1:n,1:n,sum(Tw1,2),n,n);
CO = To1 - sparse(1:n,1:n,sum(To1,2),n,n);
CWCgs = Tw1a - sparse(1:n,1:n,sum(Tw1a,2),n,n);
COCgs = To1Rs + To1a - sparse(1:n,1:n,sum(To1Rs + To1a,2),n,n);
CG = Tg1 - sparse(1:n,1:n,sum(Tg1,2),n,n);
CL = CG + CWCgs + COCgs;
