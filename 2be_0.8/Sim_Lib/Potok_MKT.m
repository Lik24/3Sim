function [TL,TW,TO,TG,TWCgs,TOCgs,TP]=Potok_MKT(T,P,KWOG,Cp,PR,rc,fp,kms,L,SGM)
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

Bw = SGM.Bwog(:,2);
Bo = SGM.Bwog(:,4);
Bg = SGM.Bwog(:,6);
Rs = SGM.Rs(:,2);

Bw1 = 0.5*(Bw(c) + Bw(r));
Bo1 = 0.5*(Bo(c) + Bo(r));
Bg1 = 0.5*(Bg(c) + Bg(r));
Rs1 = 0.5*(Rs(c) + Rs(r));

% Kfw=Sat_cal(Swe,1,1,as,aw); %water
% Kfo=Sat_cal(Swe,2,1,as,aw); %oil
% Kfp=Sat_cal(Swe.*Cpe,1,1,as,aw); %polim

[TP,wmu]=poly_vis(Cpe,PR.mup,fp,Kfw1,T,n,r,c,PR.mu(1));

% Tw=T.*Kfw1./mu(1);
Tw = T.*Kfw1./Bw1./wmu;
To = T.*Kfo1./Bo1/PR.mu(2);
Tg = T.*Kfg1./Bg1/PR.mu(3);      %//////////
 
 if kms~=0
     dPL=abs(dP./L((r+(c-1)*n)));
    %  проверить првильность задания плотности
    [Tw] = Forh(Tw,kms, PR.Ro(1), Kfw1, Ke, PR.mu(1), dPL);
    [To] = Forh(To,kms, PR.Ro(2), Kfo1, Ke, PR.mu(2), dPL);
end;

% Tw1=sparse((r+(c-1)*n),ones(size(r)),Tw,n*n,1);
% To1=sparse((r+(c-1)*n),ones(size(r)),To,n*n,1);

Twa = -Tw.*SGM.Cgsw(c);
Tw1a = sparse(r,c,Twa,n,n);     Tw1a = Tw1a + Tw1a';
Tw1 = sparse(r,c,Tw,n,n);       Tw1 = Tw1 + Tw1';

Toa = -To.*SGM.Cgso(c);
To1a = sparse(r,c,Toa,n,n);     To1a = To1a + To1a';
ToRs = To.*Rs1;
To1Rs = sparse(r,c,ToRs,n,n);   To1Rs = To1Rs + To1Rs';
To1 = sparse(r,c,To,n,n);       To1 = To1 + To1';

Tg1 = sparse(r,c,Tg,n,n);       Tg1 = Tg1 + Tg1';     %//////////

% T1=sparse(r,c,abs(dP),n,n); T1=T1+T1';
%  T1-T1'
%  Tw1-Tw1'

TW = Tw1 - sparse(1:n,1:n,sum(Tw1,2),n,n);
TO = To1 - sparse(1:n,1:n,sum(To1,2),n,n);
TWCgs = Tw1a - sparse(1:n,1:n,sum(Tw1a,2),n,n);
TOCgs = To1Rs + To1a - sparse(1:n,1:n,sum(To1Rs + To1a,2),n,n);
TG = Tg1 - sparse(1:n,1:n,sum(Tg1,2),n,n);
TL = TG + TWCgs + TOCgs;

% TP=Tp-sparse(1:n,1:n,sum(Tp,2),n,n);
% Tp1=sparse(r,c,Tp,n,n);
