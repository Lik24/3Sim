function [C,TC,WC]=GeoMex(C,TC,P,bet,P0,WC,Won)


n=size(C,1);
[r,c]=find(C);

C=C(r+(c-1)*n);

mP=(P(r)+P(c))/2;
mP0=(P0(r)+P0(c))/2;

C=re_pres1(C,mP,bet,mP0);
C=sparse(r,c,C,n,n);

U=triu(C);
[r2,c2]=find(U);
mP=(P(r2)+P(c2))/2;
mP0=(P0(r2)+P0(c2))/2;
TC=re_pres1(TC,mP,bet,mP0);
WC=re_pres1(WC,P(Won),bet,P0(Won));

end

function A=re_pres(A1,P,bet,P0)
co=(P<=P0).*(0.1)+(P>P0).*(1-bet.*(P0-P)).^3;
    A=co.*A1;
end

function A=re_pres1(A1,P,bet,P0)
%co=(P<=20).*(0.1)+(P>20).*(P<=120).*(900*P-8000)/100000+(P>120).*(1);
co=(P<=20).*(0.1)+(P>20).*(P<=60).*(2250*P-35000)/100000+(P>60).*(1);
    A=co.*A1;
end