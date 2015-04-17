function [TL,TW] = dT_dx2(T,P,Sw,Cp,as,aw,mu,rc)

n=size(P,1);

r=rc(:,1);
c=rc(:,2);

Kfw=Sat_cal(Sw,1,2,as,aw); %water
Kfo=Sat_cal(Sw,2,2,as,aw); %oil

vP=P(c,1)>P(r,1);
dP=P(c,2)-P(r,2);
% Swc=Sw(r);
% Swl=Sw(c);

Kwc=Kfw(c);
Kwl=Kfw(r);

Koc=Kfo(c);
Kol=Kfo(r);

Kfw1=Kwl.*(vP==0);
Kfo1=Kol.*(vP==0);

Twp=T.*Kfw1./mu(1).*(-dP);
Top=T.*Kfo1./mu(2).*(-dP);

Tw1=sparse(r,c,Twp,n,n);
To1=sparse(r,c,Top,n,n);
%%побочные диагонали
Twp=reshape(Tw1,n,n);
Top=reshape(To1,n,n);

%%
Kfw1=Kwc.*vP;
Kfo1=Koc.*vP;

Twp_1=T.*Kfw1./mu(1).*(dP);
Top_1=T.*Kfo1./mu(2).*(dP);

Tw1=sparse(r,c,Twp_1,n,n);
To1=sparse(r,c,Top_1,n,n);
%%главные диагонали
Twp_1=reshape(Tw1,n,n);
Top_1=reshape(To1,n,n);

TL=Top+Twp-sparse(1:n,1:n,sum(Top_1+Twp_1,2),n,n);
TW=Twp-sparse(1:n,1:n,sum(Twp_1,2),n,n);
TG=1;


end

