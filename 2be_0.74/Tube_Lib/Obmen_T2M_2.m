function [A2C]=Obmen_T2M_2(a2c,Pm,Pt,mu,MCp,CCp,r,c,KFA,KFC,AA,AC)
%KF
    n=size(Pm,1);
    m=size(Pt,1);
if isempty(r)==0

   
    kaf=KFA(:,1);
    kaof=KFA(:,2);
    kcf=KFC(:,1);
    kcof=KFC(:,2);
 
    vP=Pm(r)>=Pt(c);
%     vP
%     [Pm(r),Pt(c)]
    Kafc=kaf(r);
    Kafl=kcf(c);
    
    Kaofc=kaof(r);
    Kaofl=kcof(c);
    
    Kwe=Kafc.*vP+Kafl.*(vP==0);
    Koe=Kaofc.*vP+Kaofl.*(vP==0);

    %% 
    Cpc=MCp(r);
    Cpl=CCp(c);
    
    Cpe=Cpc.*vP+Cpl.*(vP==0);

  Afo=Koe/mu(2);  
  Afw=Kwe.*((1-Cpe)./mu(1)+Cpe/mu(4));
  ACP=Kwe.*Cpe/mu(4);
%ACW+ACO

 ACL1=Afw.*AA(r)+Afo;
 ACL2=Afw.*AC(c)+Afo;

a2cl1=a2c.*ACL1;
a2cl2=a2c.*ACL2;
a2cw=a2c.*Afw;
a2co=a2c.*Afo;
a2cp=a2c.*ACP;

A2CG=1;
    
A2C.L1=sparse(r,c,a2cl1,n,m);
A2C.L2=sparse(r,c,a2cl2,n,m);
A2C.W=sparse(r,c,a2cw,n,m);
A2C.O=sparse(r,c,a2co,n,m);
A2C.P=sparse(r,c,a2cp,n,m);

else
    
    A2C.L1=sparse(n,m);
    A2C.L2=sparse(n,m);
    A2C.W=sparse(n,m);
    A2C.P=sparse(n,m);
    A2C.G=1;
end;