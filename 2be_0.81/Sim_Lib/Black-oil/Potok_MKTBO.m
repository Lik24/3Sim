function [TW,TO,TG,TORS,TP]=Potok_MKTBO(T,P,KWOG,Cp,PR,rc,fp,kms,L,CMP)
n = size(P,1);

r = rc(:,1);
c = rc(:,2);

vP = P(c,:) - P(r,:)>0; 

Kwc = KWOG.w(c);
Kwl = KWOG.w(r);

Koc = KWOG.o(c);
Kol = KWOG.o(r);

Kgc = KWOG.g(c);     %/////////////
Kgl = KWOG.g(r);  

Cpc = Cp(c);
Cpl = Cp(r);

%Swe=Swc.*vP+Swl.*(vP==0);
Cpe = Cpc.*vP(:,1) + Cpl.*(vP(:,1)==0);
Kfw = Kwc.*vP(:,1) + Kwl.*(vP(:,1)==0);
Kfo = Koc.*vP(:,2) + Kol.*(vP(:,2)==0);
Kfg = Kgc.*vP(:,3) + Kgl.*(vP(:,3)==0);        
Rs = CMP.Rs(c,2).*vP(:,2) + CMP.Rs(r,2).*(vP(:,2)==0);    

Bw = 0.5*(CMP.Bw(c,2) + CMP.Bw(r,2));
Bo = 0.5*(CMP.Bo(c,2) + CMP.Bo(r,2));
Bg = 0.5*(CMP.Bg(c,2) + CMP.Bg(r,2));
%Rs = 0.5*(CMP.Rs(c,2) + CMP.Rs(r,2));        %///////

[TP,wmu]=poly_vis(Cpe,PR.mup,fp,Kfw,T,n,r,c,PR.mu(1));

Tw = T.*Kfw./Bw./wmu;
To = T.*Kfo./Bo/PR.mu(2);
Tg = T.*Kfg./Bg/PR.mu(3);

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
