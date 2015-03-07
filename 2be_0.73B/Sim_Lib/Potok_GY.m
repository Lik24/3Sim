function [bAl,bAw,bl,bw]=Potok_GY(T,Pgy,Pin,rc2,kfw,kfo,c,mu,na)
if na~=0
    v2=Pgy>=Pin(rc2(:,2));
    
    kfw_gy=kfw(rc2(:,1))/mu(1);
    kfo_gy=kfo(rc2(:,1))/mu(2);
    
    kfw_in=kfw(rc2(:,2))/mu(1);
    kfo_in=kfo(rc2(:,2))/mu(2);
    
    Kl=(kfw_gy+kfo_gy).*v2+(kfw_in+kfo_in).*(v2==0);
    Kw=kfw_gy.*v2+kfw_in.*(v2==0);
    
    bAl=T.*Kl;
    bAw=T.*Kw;
    
    bl=T.*Kl.*Pgy;
    bw=T.*Kw.*Pgy;
    
   
    bl = accumarray(c,bl)';
    %bl=A(v1==1)';
    bw = accumarray(c,bw)';
    %bw=A(v1==1)';
    bAl = accumarray(c,bAl)';
    %bAl=A(v1==1)';
    bAw = accumarray(c,bAw)';
    %bAw=A(v1==1)';
    
%     bl=re_bld(bl,rc2,v1,na);
%     bw=re_bld(bw,rc2,v1,na);
%     bAl=re_bld(bAl,rc2,v1,na);
%     bAw=re_bld(bAw,rc2,v1,na);
else
    bl=zeros(1,0);
    bw=zeros(1,0);
    bAl=zeros(1,0);
    bAw=zeros(1,0);
end

end

function bl=re_bld(b1,rc,v1,na)
 bl=sparse(rc(:,1),rc(:,2),b1,na,na);
 bl=sum(bl,1);
 bl=bl(v1==1);
end