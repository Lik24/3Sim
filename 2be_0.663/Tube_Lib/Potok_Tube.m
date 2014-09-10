function [CL,CW,CP,CG]=Potok_Tube(Cm,P,Sw,Cp,PR,mup,fp,kms,L)

ts=PR.ts;
tw=PR.tw;
mu=PR.mu;
Ro=PR.Ro;
Kc=PR.Kc;

n=size(Cm,1);
[r,c]=find(Cm>0);
C=Cm(r+(c-1)*n);

Kfw=Sat_tube(Sw,1,1,ts,tw); %water
Kfo=Sat_tube(Sw,2,1,ts,tw); %oil

dP=P(c)-P(r);
vP=dP>0;

Kwc=Kfw(r);
Kwl=Kfw(c);

Koc=Kfo(r);
Kol=Kfo(c);

Cpc=Cp(r);
Cpl=Cp(c);

%Swe=Swc.*vP+Swl.*(vP==0);
Cpe=Cpc.*vP+Cpl.*(vP==0);
Kfw=Kwc.*vP+Kwl.*(vP==0);
Kfo=Koc.*vP+Kol.*(vP==0);

% Kfw=Sat_tube(Swe,1,1,ts,tw); %water
% Kfo=Sat_tube(Swe,2,1,ts,tw); %oil
%Kfp=Sat_tube(Swe.*Cpe,1,1,ts,tw); %polim

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


Tw1=sparse(r,c,Tw,n,n);
To1=sparse(r,c,To,n,n);
%Tp1=sparse(r,c,Tp,n,n);

Tw=reshape(Tw1,n,n);
To=reshape(To1,n,n);
%Tp=reshape(Tp1,n,n);

CL=To+Tw-sparse(1:n,1:n,sum(To+Tw,2),n,n);
CW=Tw-sparse(1:n,1:n,sum(Tw,2),n,n);
CG=1;
