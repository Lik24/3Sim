function [A2CW,A2CO]=Obmen_T2M_2(a2c,Pm,Pt,mu,MCp,CCp,r,c,KFA,KFC,SGM)
%KF
   na = size(Pm,1);
   nc = size(Pt,1);
if isempty(r)==0

    kaf=KFA(:,1);
    kaof=KFA(:,2);
    kagf=KFA(:,3);
    
    kcf=KFC(:,1);
    kcof=KFC(:,2);
    kcgf=KFC(:,3);
 
    Bw = 0.5*(SGM.Bwog(c,2) + SGM.Bwog(r,2));
    Bo = 0.5*(SGM.Bwog(c,4) + SGM.Bwog(r,4));
    Bg = 0.5*(SGM.Bwog(c,6) + SGM.Bwog(r,6));
    Rs = 0.5*(SGM.Rs(c) + SGM.Rs(r));
    Cw = 0.5*(SGM.Cw(c) + SGM.Cw(r));
    Co = 0.5*(SGM.Co(c) + SGM.Co(r));
	Cor = 0.5*(SGM.Cor(c) + SGM.Cor(r));
	Cg = 0.5*(SGM.Cg(c) + SGM.Cg(r));
    
    vP=Pm(r)>=Pt(c);
%     vP
%     [Pm(r),Pt(c)]
    Kafc=kaf(r);
    Kafl=kcf(c);
    
    Kaofc=kaof(r);
    Kaofl=kcof(c);
    
    Kagfc=kagf(r);
    Kagfl=kcgf(c);
      
    Kwe=Kafc.*vP+Kafl.*(vP==0);
    Koe=Kaofc.*vP+Kaofl.*(vP==0);
    Kge=Kagfc.*vP+Kagfl.*(vP==0); 
    %% 
    Cpc=MCp(r);
    Cpl=CCp(c);
    
    Cpe=Cpc.*vP+Cpl.*(vP==0);
  
    ACO=Koe./Bo/mu(2);  
    ACW=Kwe.*((1-Cpe)./Bw./mu(1)+Cpe/mu(4));
    ACP=Kwe.*Cpe/mu(4);
    ACG=Kge./Bg/mu(3);
   
    A2CW = a2c.*ACW;
    A2CO = a2c.*ACO;
    A2CP = a2c.*ACP;
    A2CG = Cg.*a2c.*ACG;
    A2CWCgs = -Cw.*A2CW;
    A2COCgs = (Rs.*Cor - Co).*A2CO;
      
    A2CL = sparse(r,c,A2CW + A2CWCgs + A2COCgs,na,nc);
    A2CW = sparse(r,c,A2CW,na,nc);
    A2CO = sparse(r,c,A2CO,na,nc);
    A2CWCgs = sparse(r,c,A2CWCgs,na,nc);
    A2COCgs = sparse(r,c,A2COCgs,na,nc);
    A2CG = sparse(r,c,A2CG,na,nc);  %////////
    A2CP = sparse(r,c,A2CP,na,nc);
    
else
    A2CL = sparse(r,c,a2c,na,nc);
    A2CW = A2CL;
    A2CO = A2CL;
    A2CWCgs = A2CL;
    A2COCgs = A2CL;
    A2CG = A2CL;
end;