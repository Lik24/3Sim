function [A2CL,A2CW,A2CO,A2CP]=Obmen_T2MF2(a2c,Pm,Pt,mu,MCp,CCp,r,c,KFA,KFC,CMP,va,vc)
%KF                           
   na = size(Pm,1);
   nc = size(Pt,1);
if isempty(r)==0

    kaf=KFA(:,1);
    kaof=KFA(:,2);
    
    kcf=KFC(:,1);
    kcof=KFC(:,2);

    Bw = 0.5*(CMP.Bw(vc(c),2) + CMP.Bw(va(r),2));
    Bo = 0.5*(CMP.Bo(vc(c),2) + CMP.Bo(va(r),2));
    
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
  
    ACO=Koe./Bo/mu(2);  
    ACW=Kwe.*((1-Cpe)./Bw./mu(1)+Cpe/mu(4));
    ACP=Kwe.*Cpe/mu(4);
   
    A2CW = a2c.*ACW;
    A2CO = a2c.*ACO;
    A2CP = a2c.*ACP;
    
    A2CW = sparse(r,c,A2CW,na,nc);
    A2CO = sparse(r,c,A2CO,na,nc);
    A2CP = sparse(r,c,A2CP,na,nc);
    A2CL = sparse(1:na,1:nc,CMP.Cw(vc(c)),na,nc)*A2CW + A2CO;
    
else
    A2CW = sparse(r,c,a2c,na,nc);
    A2CO = A2CW;
    A2CL = A2CW;
end;