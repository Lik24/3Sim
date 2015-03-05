function [bAl,bAw,bAo,ql,qw,qo]=Potok_GY(T,Pgy,Pin,rc2,kfw,kfo,kfg,v1,mu,na,SGM)
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

ql = bAl.*(Pgy - Pin(rc2(:,2)));
qw = bAw.*(Pgy - Pin(rc2(:,2)));
qo = bAo.*(Pgy - Pin(rc2(:,2)));

ql=re_bld(ql,rc2,v1,na);
qw=re_bld(qw,rc2,v1,na);
qo=re_bld(qo,rc2,v1,na);
bAl=re_bld(bAl,rc2,v1,na);
bAw=re_bld(bAw,rc2,v1,na);
bAo=re_bld(bAo,rc2,v1,na);
ql = ql';
qw = qw';
qo = qo';
end

function bl=re_bld(b1,rc,v1,na)
 bl=sparse(rc(:,1),rc(:,2),b1,na,na);
 bl=sum(bl,1);
 bl=bl(v1==1);
end