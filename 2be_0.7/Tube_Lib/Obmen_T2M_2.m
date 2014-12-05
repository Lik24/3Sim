function [A2CL,A2CW,A2CP,A2CG]=Obmen_T2M_2(a2c,Pm,Pt,mu,MCp,CCp,r,c,KFA,KFC)
%KF

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

  ACO=Koe/mu(2);  
  ACW=Kwe.*((1-Cpe)./mu(1)+Cpe/mu(4));
  ACP=Kwe.*Cpe/mu(4);
%ACW+ACO

    A2CL=a2c.*(ACW+ACO);
    A2CW=a2c.*ACW;
    A2CP=a2c.*ACP;
    A2CG=1;
else
    
    A2CL=a2c;
    A2CW=a2c;
    A2CP=a2c;
    
    A2CG=1;
end;