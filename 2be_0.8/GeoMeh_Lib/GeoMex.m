function [C,TC]=GeoMex(C,TC,P,fll)

%C
n=size(C,1);
[r,c]=find(C);
[r1,c1]=find(C>0);

C=C(r+(c-1)*n);
mP=(P(r)+P(c))/2;

C=re_pres(C,mP,fll);
C=sparse(r,c,C,n,n);

mP=(P(r1)+P(c1))/2;
TC=re_pres(TC,mP,fll);

end

function A=re_pres(A1,P,fll)
if fll==10
    co=ones(size(P));
elseif fll<5
    co=110*log(P)-429.3;
elseif fll==3
    co=0.01*P-0.5;
elseif fll==4
     co=0.01*P-0.5;
elseif fll==5
     co=0.01*P-0.5;
elseif fll==6
     co=0.01*P-0.5;
elseif fll==7
     co=0.01*P-0.5;
end;
% size(co)
% size(A1)
    A=co.*A1;
end
