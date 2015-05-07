function [bAl,bAw,bAo,ql,qw,qo]=Potok_GY(T,Pgy,Pin,dPin,rc2,kfw,kfo,kfg,v1,mu,na,c,SGM,ql,qw,qo)
v2=Pgy>=Pin(rc2(:,2));

Bw = 0.5*(SGM.Bwog(rc2(:,2),2) + SGM.Bwog(rc2(:,1),2));
Bo = 0.5*(SGM.Bwog(rc2(:,2),4) + SGM.Bwog(rc2(:,1),4));
Bg = 0.5*(SGM.Bwog(rc2(:,2),6) + SGM.Bwog(rc2(:,1),6));
Rs = 0.5*(SGM.Rs(rc2(:,2)) + SGM.Rs(rc2(:,1)));

Kw = (kfw(rc2(:,1)).*v2+kfw(rc2(:,2)).*(v2==0))/mu(1)./Bw;
Ko = (kfo(rc2(:,1)).*v2+kfo(rc2(:,2)).*(v2==0))/mu(2)./Bo;
Kg = (kfg(rc2(:,1)).*v2+kfg(rc2(:,2)).*(v2==0))/mu(3)./Bg;

bAw = T.*Kw;
bAo = T.*Ko;
bAl = T.*Kg - SGM.Cgsw(rc2(:,2)).*bAw + (Rs - SGM.Cgso(rc2(:,2))).*bAo;

ql1 = - bAl.*dPin(rc2(:,2));
qw1 = - bAw.*dPin(rc2(:,2));
qo1 = - bAo.*dPin(rc2(:,2));

ql1=re_bld(ql1,rc2,v1,na);
ql = ql + ql1';
qw1=re_bld(qw1,rc2,v1,na);
qw = qw + qw1';
qo1=re_bld(qo1,rc2,v1,na);
qo = qo + qo1';

%bAl=re_bld(bAl,rc2,v1,na);
%bAw=re_bld(bAw,rc2,v1,na);
%bAo=re_bld(bAo,rc2,v1,na);

bAl = accumarray(c,bAl)';
bAw = accumarray(c,bAw)';
bAo = accumarray(c,bAo)';
end

function bl=re_bld(b1,rc,v1,na)
 bl=sparse(rc(:,1),rc(:,2),b1,na,na);
 bl=sum(bl,1);
 bl=bl(v1==1);
end