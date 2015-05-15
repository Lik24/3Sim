function [TW,TO,TG,TORS,Twl,Tol,Trsl,Tgl,TSw,TSo,TPb,ToSw,TRoSw,TP]=Potok_MKTSeqBO(T,P,KWOG,Cp,PR,rc,fp,kms,L,CMP,SGM)
n = size(P,1);

r = rc(:,1);
c = rc(:,2);

vP = P(c,:) - P(r,:)>0; 

Cpe = Cp(c).*vP(:,1) + Cp(r).*(vP(:,1)==0);
Kfw = KWOG.w(c).*vP(:,1) + KWOG.w(r).*(vP(:,1)==0); dKfw = KWOG.dwds(c).*vP(:,1) + KWOG.dwds(r).*(vP(:,1)==0);
Kfo = KWOG.o(c).*vP(:,2) + KWOG.o(r).*(vP(:,2)==0); dKfow = KWOG.dodsw(c).*vP(:,2) + KWOG.dodsw(r).*(vP(:,2)==0);
Kfg = KWOG.g(c).*vP(:,3) + KWOG.g(r).*(vP(:,3)==0); dKfog = KWOG.dodsg(c).*vP(:,2) + KWOG.dodsg(r).*(vP(:,2)==0);  
dKfg = KWOG.dgds(c).*vP(:,3) + KWOG.dgds(r).*(vP(:,3)==0); 
Rs = CMP.Rs(c,2).*vP(:,2) + CMP.Rs(r,2).*(vP(:,2)==0);    

Bw = 0.5*(CMP.Bw(c,2) + CMP.Bw(r,2));
Bo = 0.5*(CMP.Bo(c,2) + CMP.Bo(r,2));
Bg = 0.5*(CMP.Bg(c,2) + CMP.Bg(r,2));
dBwP = 0.5*(SGM.dBwP(c) + SGM.dBwP(r));
dBoP = 0.5*(SGM.dBoP(c) + SGM.dBoP(r));
dBoPb = 0.5*(SGM.dBoPb(c) + SGM.dBoPb(r));
dBgP = 0.5*(SGM.dBgP(c) + SGM.dBgP(r));

[TP,wmu]=poly_vis(Cpe,PR.mup,fp,Kfw,T,n,r,c,PR.mu(1));

Tw = T./Bw./wmu;
To = T./Bo/PR.mu(2);
Tg = T./Bg/PR.mu(3);

Tws = Tw.*dKfw.*(P(c,1)-P(r,1));
Toso = -To.*dKfog.*(P(c,2)-P(r,2));
Tosw = To.*(dKfow-dKfog).*(P(c,2)-P(r,2));

Tw = Tw.*Kfw; 
To = To.*Kfo;
Tg = Tg.*Kfg;

Tpb = To.*(PR.rs - Rs.*dBoPb./Bo).*(P(c,2)-P(r,2));

 if kms~=0
     dPL=abs(dP./L((r+(c-1)*n)));
    %  проверить првильность задания плотности
    [Tw] = Forh(Tw,kms, PR.Ro(1), Kfw1, Ke, PR.mu(1), dPL);
    [To] = Forh(To,kms, PR.Ro(2), Kfo1, Ke, PR.mu(2), dPL);
end;

% Tw1=sparse((r+(c-1)*n),ones(size(r)),Tw,n*n,1);
% To1=sparse((r+(c-1)*n),ones(size(r)),To,n*n,1);

Tw1 = sparse(r,c,-0.5*Tw.*dBwP./Bw.*(P(c,1)-P(r,1)),n,n);    Tw1 = Tw1 - Tw1';
To1 = sparse(r,c,-0.5*To.*dBoP./Bo.*(P(c,2)-P(r,2)),n,n);       To1 = To1 - To1';
ToRs1 = sparse(r,c,0.5*To.*(PR.rs-Rs.*dBoP./Bo).*(P(c,2)-P(r,2)),n,n); ToRs1 = ToRs1 - ToRs1';
Tg1 = sparse(r,c,-0.5*Tg.*dBgP./Bg.*(P(c,3)-P(r,3)),n,n);       Tg1 = Tg1 - Tg1';  
Twl = Tw1 - sparse(1:n,1:n,sum(Tw1,2),n,n);
Tol = To1 - sparse(1:n,1:n,sum(To1,2),n,n);
Trsl = ToRs1 - sparse(1:n,1:n,sum(ToRs1,2),n,n);
Tgl = Tg1 - sparse(1:n,1:n,sum(Tg1,2),n,n);

TSw = sparse(r,c,Tws.*vP(:,1),n,n) + sparse(c,r,-Tws.*(~vP(:,1)),n,n);    
TSo = sparse(r,c,Toso.*vP(:,2),n,n) + sparse(c,r,-Toso.*(~vP(:,2)),n,n);
TPb = sparse(r,c,Tpb.*vP(:,2),n,n) + sparse(c,r,-Tpb.*(~vP(:,2)),n,n);
TSw = TSw - sparse(1:n,1:n,sum(TSw,1),n,n);
TPb = TPb - sparse(1:n,1:n,sum(TPb,1),n,n);
TSo = TSo - sparse(1:n,1:n,sum(TSo,1),n,n);

ToSw = sparse(r,c,Tosw.*vP(:,2),n,n) + sparse(c,r,-Tosw.*(~vP(:,2)),n,n);
TRoSw = sparse(r,c,Tosw.*Rs.*vP(:,2),n,n) + sparse(c,r,-Tosw.*Rs.*(~vP(:,2)),n,n);
ToSw = ToSw - sparse(1:n,1:n,sum(ToSw,1),n,n);
TRoSw = TRoSw - sparse(1:n,1:n,sum(TRoSw,1),n,n);

Tw = sparse(r,c,Tw,n,n);    Tw = Tw + Tw';
ToRs = sparse(r,c,To.*Rs,n,n);     ToRs = ToRs + ToRs';
To = sparse(r,c,To,n,n);       To = To + To';
Tg = sparse(r,c,Tg,n,n);       Tg = Tg + Tg';  

TW = Tw - sparse(1:n,1:n,sum(Tw,2),n,n);
TO = To - sparse(1:n,1:n,sum(To,2),n,n);
TORS = ToRs - sparse(1:n,1:n,sum(ToRs,2),n,n);
TG = Tg - sparse(1:n,1:n,sum(Tg,2),n,n);


% TP=Tp-sparse(1:n,1:n,sum(Tp,2),n,n);
% Tp1=sparse(r,c,Tp,n,n);
