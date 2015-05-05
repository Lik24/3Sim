function [A2CL1,A2CL2,A2CW,A2CP,A2CG]=Obmen_T2M(A2C,Pm,Pt,A_kfw,A_kfo,C_kfw,C_kfo,K,PR,MCp,CCp,AA,AC)

if isempty(Pm)==0 && isempty(Pt)==0
mu=PR.mu;
mup=PR.mup;

[n,m]=size(A2C);
[r,c]=find(A2C);
A2C=A2C(:);
a2c=A2C(r+(c-1)*n).*K(r);

% Pm=P(1:n);
% Pt=P(n+1:end);

vP=Pm(r)>=Pt(c);
Cpc=MCp(r);
Cpl=CCp(c);

Cpe=Cpc.*vP+Cpl.*(vP==0);
kfw=A_kfw(r).*vP+C_kfw(c).*(vP==0);
kfo=A_kfo(r).*vP+C_kfo(c).*(vP==0);

wmu=fusion(Cpe,mup);

Afw=kfw./wmu; 
Afp=kfw./wmu.*Cpe; 
Afo=kfo/mu(2); %oil

ACL1=Afw.*AA(r)+Afo;
ACL2=Afw.*AC(c)+Afo;
ACW=Afw;
ACP=Afp;
%ACL
a2cl1=a2c.*ACL1;
a2cl2=a2c.*ACL2;
a2cw=a2c.*ACW;
a2cp=a2c.*ACP;

A2CL1=sparse(r,c,a2cl1,n,m);
A2CL2=sparse(r,c,a2cl2,n,m);
A2CW=sparse(r,c,a2cw,n,m);
A2CP=sparse(r,c,a2cp,n,m);
A2CG=1;
else
    A2CL1=A2C;
    A2CL2=A2C;
    A2CW=A2C;
    A2CP=A2C;
    A2CG=1;  
end