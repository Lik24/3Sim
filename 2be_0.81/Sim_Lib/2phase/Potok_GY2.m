function [bAl,ql,qw]=Potok_GY2(T,Pgy,Pin,dPin,rc2,kfw,kfo,mu,c,CMP,qw,ql)
v2=Pgy>=Pin(rc2(:,2));

Bw = 0.5*(CMP.Bw(rc2(:,2),2) + CMP.Bw(rc2(:,1),2));
Bo = 0.5*(CMP.Bo(rc2(:,2),2) + CMP.Bo(rc2(:,1),2));

Kw = (kfw(rc2(:,1)).*v2+kfw(rc2(:,2)).*(v2==0))/mu(1)./Bw;
Ko = (kfo(rc2(:,1)).*v2+kfo(rc2(:,2)).*(v2==0))/mu(2)./Bo;

bAw = T.*Kw;
bAo = T.*Ko;
bAl = CMP.Cw(rc2(:,2)).*bAw + bAo;

ql1 = - bAl.*dPin(rc2(:,2));
qo1 = - bAw.*dPin(rc2(:,2));
qw1 = - bAw.*dPin(rc2(:,2));

qw = qw + accumarray(c,qw1);
ql = ql + accumarray(c,ql1);
bAl = accumarray(c,bAl);
end

function bl=re_bld(b1,rc,v1,na)
 bl=sparse(rc(:,1),rc(:,2),b1,na,na);
 bl=sum(bl,1);
 bl=bl(v1==1);
end