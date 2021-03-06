function [dF1,dF2]=dif_F_dx2(T,P,Cp,dKfw,dKfo,mu,rc,Uf,CpW,Pw,WonX,D)

n=size(P,1);

Won=WonX(:,1);
Wf=WonX(:,2);
Uf=Uf(WonX(:,3));
CpW=CpW(WonX(:,3));

r=rc(:,1);
c=rc(:,2);

vP=P(c,1)>P(r,1);
dP=P(c,1)-P(r,1);
% Swc=Sw(r);
% Swl=Sw(c);

Kwc=dKfw(c);
Kwl=dKfw(r);

Koc=dKfo(c);
Kol=dKfo(r);

Kfw1=Kwl.*(vP==0);
Kfo1=Kol.*(vP==0);

Twp=T.*Kfw1./mu(1).*(-dP);
Top=T.*Kfo1./mu(2).*(-dP);

Twp=sparse(r,c,Twp,n,n);
Top=sparse(r,c,Top,n,n);

%%
Kfw1=Kwc.*vP;
Kfo1=Koc.*vP;

Twp_1=T.*Kfw1./mu(1).*(dP);
Top_1=T.*Kfo1./mu(2).*(dP);

Tw1=sparse(r,c,Twp_1,n,n);
To1=sparse(r,c,Top_1,n,n);

dTL=Top+Twp-sparse(1:n,1:n,sum(To1+Tw1,2),n,n);
dTW=Twp-sparse(1:n,1:n,sum(Tw1,2),n,n);


[dW1,dW6,~]=Well_MKT_2(Wf,Won,Uf,Cp,mu,CpW,dKfw,dKfo);  %производные
dPw=Pw(WonX(:,3))-P(Won);

dF1=dTL+sparse(Won,Won,dW1.*dPw,n,n);
dF2=D+dTW+sparse(Won,Won,dW6.*dPw,n,n);

