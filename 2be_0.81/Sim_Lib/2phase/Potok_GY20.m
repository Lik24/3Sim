function [bAl,ql,qo,qw]=Potok_GY20(T,Pgy,Pin,rc2,kfw,kfo,v1,mu,na,c,CMP)
v2=Pgy>=Pin(rc2(:,2));

Bw = 0.5*(CMP.Bw(rc2(:,2),2) + CMP.Bw(rc2(:,1),2));
Bo = 0.5*(CMP.Bo(rc2(:,2),2) + CMP.Bo(rc2(:,1),2));

Kw = (kfw(rc2(:,1)).*v2+kfw(rc2(:,2)).*(v2==0))/mu(1)./Bw;
Ko = (kfo(rc2(:,1)).*v2+kfo(rc2(:,2)).*(v2==0))/mu(2)./Bo;

bAw = T.*Kw;
bAo = T.*Ko;
bAl = CMP.Cw(rc2(:,2)).*bAw + bAo;

ql = bAl.*(Pgy - Pin(rc2(:,2)));
qo = bAo.*(Pgy - Pin(rc2(:,2)));
qw = bAw.*(Pgy - Pin(rc2(:,2)));

bAl = accumarray(c,bAl);
ql = accumarray(c,ql);
qo = accumarray(c,qo);
qw = accumarray(c,qw);
end

function bl=re_bld(b1,rc,v1,na)
 bl=sparse(rc(:,1),rc(:,2),b1,na,na);
 bl=sum(bl,1);
 bl=bl(v1==1);
end