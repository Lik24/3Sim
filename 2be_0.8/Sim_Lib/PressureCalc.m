function [Pi,CL,GL,QWOG,Phi,SGMMD,Qm2,Qd2]=PressureCalc(Pi,Sw,So,Phi,Sw0,So0,Cp,TRM,KWOG,SGM,RC,WELL,fp,VEC,GEOM,DATA,NCELL,M2FR,ft,PR,QWOG,BOUNDFL,dt,t,Qf,QBOUND,Qz)

  na = NCELL.na;
  nc = NCELL.nc;
  ng = NCELL.ng;
  nd = NCELL.nd;
  nb = NCELL.nb;
  nw = NCELL.nw;
  Uf = WELL.Uf(:,ft+1);
  CpW = WELL.CpW(:,ft+1);
  b1wb=sparse(nb,1);
  vad = RC.ADr;
  dPw = zeros(size(Qz(:,1),1),1);
  Qm2=zeros(nw,5);
  Qd2=zeros(nw,5);
  kj = 0;
  flag_gim = 1;
  while flag_gim==1 && kj<10
    kj=kj+1;
   
   [TL,TW,TO,TGG,TWCgs,TOCgs,TP]=Potok_MKT(TRM.TTM,Phi(VEC.va,:),KWOG,Cp(VEC.va,1),PR,RC.Arc2,fp,PR.kms(1),GEOM.L,SGM);  %проводимости по фазам
   [CL,CW,CO,CG,CWCgs,COCgs,~]=Potok_Tube(TRM.TC,Phi(VEC.vc,:),Sw(VEC.vc,t),Cp(VEC.vc,t),PR,fp,PR.kms(2),DATA.Lc,RC.Cr2,RC.Cc2,NCELL.nc,SGM);
   [GL,GW,GO,GG,GWCgs,GOCgs,~]=Potok_Tube(TRM.TG,Phi(VEC.vg,:),Sw(VEC.vg,t),Cp(VEC.vg,t),PR,fp,PR.kms(3),GEOM.Lg,RC.Gr2,RC.Gc2,NCELL.ng,SGM);
   [DL,DW,DO,DG,DWCgs,DOCgs,DP]=Potok_Tube(TRM.TD,Phi(VEC.vd,:),Sw(VEC.vd,1),Cp(VEC.vd,1),PR,fp,PR.kms(4),GEOM.Ld,RC.Dr2,RC.Dc2,NCELL.nd,SGM);
                                                   
   [A2CL,A2CW,A2CO,A2CG,A2CWCgs,A2COCgs,~]=Obmen_T2M(M2FR.A2C,Pi(VEC.va,1),Pi(VEC.vc,1),Sw(VEC.va,1),Sw(VEC.vc,t),KWOG.K(:,1),PR,Cp(VEC.va,1),Cp(VEC.vc,t),SGM);
   [A2GL,A2GW,A2GO,A2GG,A2GWCgs,A2GOCgs,~]=Obmen_T2M(M2FR.A2G,Pi(VEC.va,1),Pi(VEC.vg,1),Sw(VEC.va,1),Sw(VEC.vg,t),KWOG.K(:,1),PR,Cp(VEC.va,1),Cp(VEC.vg,t),SGM);
        
   [A2DL,A2DW,A2DO,A2DG,A2DWCgs,A2DOCgs,A2DP]=Obmen_T2M(M2FR.A2D,Pi(VEC.va,1),Pi(VEC.vd,1),Sw(VEC.va,1),Sw(VEC.vd,1),KWOG.K(:,1),PR,Cp(VEC.va,1),Cp(VEC.vd,1),SGM);
   [A2BL,A2BW,A2BO,A2BG,A2BWCgs,A2BOCgs,A2BP]=Obmen_T2M(M2FR.A2B,Pi(VEC.va,1),Pi(VEC.vb,1),Sw(VEC.va,1),Sw(VEC.vb,1),ones(na,1),PR,Cp(VEC.va,1),Cp(VEC.vb,1),SGM);
   [D2BL,D2BW,D2BO,D2BG,D2BWCgs,D2BOCgs,D2BP]=Obmen_T2M(M2FR.D2B,Pi(VEC.vd,1),Pi(VEC.vb,1),Sw(VEC.vd,1),Sw(VEC.vb,1),KWOG.K(:,1),PR,Cp(VEC.vd,1),Cp(VEC.vb,1),SGM);
        
   [D2CL,D2CW,D2CO,D2CG,D2CWCgs,D2COCgs,~]=Obmen_T2M(M2FR.C2D,Pi(VEC.vd,1),Pi(VEC.vc,1),Sw(VEC.vd,1),Sw(VEC.vc,t),KWOG.K(:,1),PR,Cp(VEC.vd,1),Cp(VEC.vc,t),SGM);
   [D2GL,D2GW,D2GO,D2GG,D2GWCgs,D2GOCgs,~]=Obmen_T2M(M2FR.G2D,Pi(VEC.vd,1),Pi(VEC.vg,1),Sw(VEC.vd,1),Sw(VEC.vg,t),KWOG.K(:,1),PR,Cp(VEC.vd,1),Cp(VEC.vg,t),SGM);
        
        %Wf=Wf0.*(1-0.0001*(P0(Won)-Pi(Won))).^3;
   [W1,W6,Wo,Wg,W7]=Well_MKT(WELL.Won,Uf(WELL.Won(:,3)),Sw(VEC.va,1),Cp(VEC.va,1),PR.aw,PR.as,PR,CpW(WELL.Won(:,3)),SGM);% W1 - проводимость по всей жидкости, W6 - только для воды, W7 - полимер 
   [W1C,W6C,WoC,WgC,W7C]=Well_MKT(WELL.WonC,Uf(WELL.WonC(:,3)),Sw(VEC.vc,t),Cp(VEC.vc,t),PR.tw,PR.ts,PR,CpW(WELL.WonC(:,3)),SGM);
   [W1G,W6G,WoG,WgG,W7G]=Well_MKT(WELL.WonG,Uf(WELL.WonG(:,3)),Sw(VEC.vg,t),Cp(VEC.vg,t),PR.tw,PR.ts,PR,CpW(WELL.WonG(:,3)),SGM);
   [W1D,W6D,WoD,WgD,W7D]=Well_MKT(WELL.WonD,Uf(WELL.WonD(:,3)),Sw(VEC.vd,1),Cp(VEC.vd,1),PR.tw,PR.ts,PR,CpW(WELL.WonD(:,3)),SGM);
        
   A1= -sparse(WELL.Won(:,1),WELL.Won(:,1),W1,na,na)-sparse(1:na,1:na,SGM.Clp(VEC.va)+sum(BOUNDFL.b1gm(:,1:2),2),na,na);  %Матрица коэф. для пор
   C1= -sparse(1:nc,1:nc,SGM.Clp(VEC.vc)',nc,nc)-sparse(WELL.WonC(:,1),WELL.WonC(:,1),W1C,nc,nc);                       %Матрица коэф. для вертикальных трещ.
   G1= -sparse(1:ng,1:ng,SGM.Clp(VEC.vg)',ng,ng)-sparse(WELL.WonG(:,1),WELL.WonG(:,1),W1G,ng,ng);                       %Матрица коэф. для гориз. трещ.
   D1= -sparse(WELL.WonD(:,1),WELL.WonD(:,1),W1D,nd,nd)-sparse(1:nd,1:nd,SGM.Clp(VEC.vd)+sum(BOUNDFL.b1gd(:,1:2),2),nd,nd);                       %Матрица коэф. для двойной пор.
   B1= -sparse(1:nb,1:nb,sum(A2BL,1)+sum(D2BL,1)+SGM.Clp(VEC.vb)'+BOUNDFL.b1gb',nb,nb);  %Матрица коэф. для границ
        
   AG1=TGG-sparse(1:na,1:na,sum(A2CG,2)+sum(A2GG,2)+sum(A2BG,2)+sum(A2DG,2),na,na);  
   CG1=CG-sparse(1:nc,1:nc,sum(A2CG,1)+sum(D2CG,1),nc,nc);                       
   GG1=GG-sparse(1:ng,1:ng,sum(A2GG,1)+sum(D2GG,1),ng,ng);   
   DG1=DG-sparse(1:nd,1:nd,sum(A2DG,1)'+sum(D2BG,2)+sum(D2CG,2)+sum(D2GG,2),nd,nd);
        
   AW1Cgs=TWCgs-sparse(1:na,1:na,sum(A2CWCgs,2)+sum(A2GWCgs,2)+sum(A2BWCgs,2)+sum(A2DWCgs,2),na,na);  
   CW1Cgs=CWCgs-sparse(1:nc,1:nc,sum(A2CWCgs,1)+sum(D2CWCgs,1),nc,nc);                       
   GW1Cgs=GWCgs-sparse(1:ng,1:ng,sum(A2GWCgs,1)+sum(D2GWCgs,1),ng,ng);   
   DW1Cgs=DWCgs-sparse(1:nd,1:nd,sum(A2DWCgs,1)'+sum(D2BWCgs,2)+sum(D2CWCgs,2)+sum(D2GWCgs,2),nd,nd);
      
   AO1Cgs=TOCgs-sparse(1:na,1:na,sum(A2COCgs,2)+sum(A2GOCgs,2)+sum(A2BOCgs,2)+sum(A2DOCgs,2),na,na);  
   CO1Cgs=COCgs-sparse(1:nc,1:nc,sum(A2COCgs,1)+sum(D2COCgs,1),nc,nc);                       
   GO1Cgs=GOCgs-sparse(1:ng,1:ng,sum(A2GOCgs,1)+sum(D2GOCgs,1),ng,ng);   
   DO1Cgs=DOCgs-sparse(1:nd,1:nd,sum(A2DOCgs,1)'+sum(D2BOCgs,2)+sum(D2COCgs,2)+sum(D2GOCgs,2),nd,nd);
   
   AW1=TW-sparse(1:na,1:na,sum(A2CW,2)+sum(A2GW,2)+sum(A2BW,2)+sum(A2DW,2),na,na);  
   CW1=CW-sparse(1:nc,1:nc,sum(A2CW,1)+sum(D2CW,1),nc,nc);                       
   GW1=GW-sparse(1:ng,1:ng,sum(A2GW,1)+sum(D2GW,1),ng,ng);   
   DW1=DW-sparse(1:nd,1:nd,sum(A2DW,1)'+sum(D2BW,2)+sum(D2CW,2)+sum(D2GW,2),nd,nd);
      
   AO1=TO-sparse(1:na,1:na,sum(A2CO,2)+sum(A2GO,2)+sum(A2BO,2)+sum(A2DO,2),na,na);  
   CO1=CO-sparse(1:nc,1:nc,sum(A2CO,1)+sum(D2CO,1),nc,nc);                       
   GO1=GO-sparse(1:ng,1:ng,sum(A2GO,1)+sum(D2GO,1),ng,ng);   
   DO1=DO-sparse(1:nd,1:nd,sum(A2DO,1)'+sum(D2BO,2)+sum(D2CO,2)+sum(D2GO,2),nd,nd);
   
   AMG=[AG1,   A2CG, A2GG, A2DG, A2BG;
       A2CG',  CG1,  M2FR.C2G, D2CG', M2FR.C2B;
       A2GG', M2FR.C2G', GG1,  D2GG', M2FR.G2B;
       A2DG', D2CG, D2GL,  DG1,  D2BG;
       A2BG', M2FR.C2B',M2FR.G2B',D2BG', zeros(nb,nb)];
   
   AMWCgs=[AW1Cgs,   A2CWCgs, A2GWCgs, A2DWCgs, A2BWCgs;
       A2CWCgs',  CW1Cgs,  M2FR.C2G, D2CWCgs', M2FR.C2B;
       A2GWCgs', M2FR.C2G', GW1Cgs,  D2GWCgs', M2FR.G2B;
       A2DWCgs', D2CWCgs, D2GWCgs,  DW1Cgs,  D2BWCgs;
       A2BWCgs', M2FR.C2B',M2FR.G2B',D2BWCgs', B1];
   
   AMOCgs=[AO1Cgs,   A2COCgs, A2GOCgs, A2DOCgs, A2BOCgs;
       A2COCgs',  CO1Cgs,  M2FR.C2G, D2COCgs', M2FR.C2B;
       A2GOCgs', M2FR.C2G', GO1Cgs,  D2GOCgs', M2FR.G2B;
       A2DOCgs', D2COCgs, D2GOCgs,  DO1Cgs,  D2BOCgs;
       A2BOCgs', M2FR.C2B',M2FR.G2B',D2BOCgs', zeros(nb,nb)];
   
   AMW=[AW1,   A2CW, A2GW, A2DW, A2BW;
       A2CW',  CW1,  M2FR.C2G, D2CW', M2FR.C2B;
       A2GW', M2FR.C2G', GW1,  D2GW', M2FR.G2B;
       A2DW', D2CW, D2GW,  DW1,  D2BW;
       A2BW', M2FR.C2B',M2FR.G2B',D2BW', B1];
   
   AMO=[AO1,   A2CO, A2GO, A2DO, A2BO;
       A2CO',  CO1,  M2FR.C2G, D2CO', M2FR.C2B;
       A2GO', M2FR.C2G', GO1,  D2GO', M2FR.G2B;
       A2DO', D2CO, D2GO,  DO1,  D2BO;
       A2BO', M2FR.C2B',M2FR.G2B',D2BO', zeros(nb,nb)];

   W2M=sparse(WELL.Won(:,3),WELL.Won(:,1),W1,nw,na);
   W2C=sparse(WELL.WonC(:,3),WELL.WonC(:,1),W1C,nw,nc);
   W2G=sparse(WELL.WonG(:,3),WELL.WonG(:,1),W1G,nw,ng);
   W2D=sparse(WELL.WonD(:,3),WELL.WonD(:,1),W1D,nw,nd);
   W2B=sparse(nw,nb);

   WM1=[W2M,W2C,W2G,W2D,W2B];
   WM2=WM1';
   W3vec=sparse(WELL.Won(:,3),1,W1,nw,1)+sparse(WELL.WonC(:,3),1,W1C,nw,1)+sparse(WELL.WonG(:,3),1,W1G,nw,1)+sparse(WELL.WonD(:,3),1,W1D,nw,1);
   WM3=-sparse(1:nw,1:nw,W3vec,nw,nw);
           
 %         while  flag_pwq==1
     %       qfSize = size(Qf~=0);
     %       PwNl=repmat(Pw(:,ft+1),Nl,1);
            % PwNl=PwNl(ka1==1);
          
   [dItime,dItimeW,dItimeO] = NR_Time(GEOM.dV,Sw,So,Sw0,So0,SGM,VEC.va,VEC.vc,VEC.vg,VEC.vd,VEC.vb,dt);         
   dPhi = AMWCgs*Phi(:,1) + AMOCgs*Phi(:,2) + AMG*Phi(:,3);

 %  while  flag_pwq==1

     b1wm=sparse(WELL.Won(:,1),ones(1,size(WELL.Won,1)),-W1.*dPw(WELL.Won(:,3)),na,1);
     b1wc=sparse(WELL.WonC(:,1),ones(1,size(WELL.WonC,1)),-W1C.*dPw(WELL.WonC(:,3),ft+1),nc,1);
     b1wg=sparse(WELL.WonG(:,1),ones(1,size(WELL.WonG,1)),-W1G.*dPw(WELL.WonG(:,3),ft+1),ng,1);
     b1wd=sparse(WELL.WonD(:,1),ones(1,size(WELL.WonD,1)),-W1D.*dPw(WELL.WonD(:,3)),nd,1);
            
    b1wm=b1wm.*(sum(W2M(Qf==0,:),1)~=0)';
    b1wc=b1wc.*(sum(W2C(Qf==0,:),1)~=0)';
    b1wg=b1wg.*(sum(W2G(Qf==0,:),1)~=0)';
    b1wd=b1wd.*(sum(W2D(Qf==0,:),1)~=0)';
            
    b1wm = b1wm - sparse(WELL.Won(:,1),ones(1,size(WELL.Won,1)),QWOG.QQm,na,1);
    b1wc = b1wc - sparse(WELL.WonC(:,1),ones(1,size(WELL.WonC,1)),QWOG.QQc,nc,1);
    b1wg = b1wg - sparse(WELL.WonG(:,1),ones(1,size(WELL.WonG,1)),QWOG.QQg,ng,1);
    b1wd = b1wd - sparse(WELL.WonD(:,1),ones(1,size(WELL.WonD,1)),QWOG.QQd,nd,1);
                           
    BM=[b1wm;b1wc;b1wg;b1wd;b1wb] - [QBOUND.Qm;QBOUND.Qc;QBOUND.Qg;QBOUND.Qd;QBOUND.Qb] + dItime - dPhi;
  
    W2M1=WM1(Qf~=0,:);
    W2M2=WM2(:,Qf~=0);
    W2M3=WM3(Qf~=0,Qf~=0);
                  
    AM = [A1, zeros(na,nc), zeros(na,ng), zeros(na,nd), zeros(na,nb);
          zeros(nc,na), C1, zeros(nc,ng), zeros(nc,nd), zeros(nc,nb);
          zeros(ng,na), zeros(ng,nc), G1, zeros(ng,nd), zeros(ng,nb);
          zeros(nd,na), zeros(nd,nc), zeros(nd,ng), D1, zeros(nd,nb);
          zeros(nb,na), zeros(nb,nc), zeros(nb,ng), zeros(nb,nd), B1];

    AM = AM + AMWCgs + AMOCgs + AMG; 
    
    Qz1 = sparse(WELL.Won(Qf~=0,1),ones(1,size(WELL.Won(Qf~=0),1)),Qz(Qf~=0,ft+1),nw,1);
    Qm1 = sparse(WELL.Won(Qf~=0,1),ones(1,size(WELL.Won(Qf~=0),1)),QWOG.QQm(Qf~=0),nw,1);
    Qc1 = sparse(WELL.WonC(Qf(WELL.WonC(:,3))~=0,3),ones(1,size(WELL.WonC(Qf(WELL.WonC(:,3))~=0,3),1)),QWOG.QQc(Qf(WELL.WonC(:,3))~=0),nw,1);
    Qg1 = sparse(WELL.WonG(Qf(WELL.WonG(:,3))~=0,3),ones(1,size(WELL.WonG(Qf(WELL.WonG(:,3))~=0,3),1)),QWOG.QQg(Qf(WELL.WonG(:,3))~=0),nw,1);
    Qd1 = sparse(WELL.WonD(Qf(WELL.WonD(:,3))~=0,3),ones(1,size(WELL.WonD(Qf(WELL.WonD(:,3))~=0,3),1)),QWOG.QQd(Qf(WELL.WonD(:,3))~=0),nw,1);
   
    Q1 = -Qz1+Qm1+Qc1+Qg1+Qd1;
    
    dPt = [BM',Q1(Qf~=0)']/[AM,W2M2;W2M1,W2M3];  %/////////////
    dPt = dPt';
   
    dPw(Qf~=0) = dPt(na+nc+ng+nd+nb+1:end);
   
    [QWOG.QQm,QQmwo]=QIter(QWOG.QQm,W6,Wo,Wg,dPt(VEC.va),WELL.Won(:,1),dPw(WELL.Won(:,3)),SGM);
    [QWOG.QQc,QQcwo]=QIter(QWOG.QQc,W6C,WoC,WgC,dPt(VEC.vc),WELL.WonC(:,1),dPw(WELL.WonC(:,3)),SGM);
    [QWOG.QQg,QQgwo]=QIter(QWOG.QQg,W6G,WoG,WgG,dPt(VEC.vg),WELL.WonG(:,1),dPw(WELL.WonG(:,3)),SGM);
    [QWOG.QQd,QQdwo]=QIter(QWOG.QQd,W6D,WoD,WgD,dPt(VEC.vd),WELL.WonD(:,1),dPw(WELL.WonD(:,3)),SGM); 
   
    QWOG.QQmwo(:,1) = QWOG.QQmwo(:,1) + sparse(WELL.Won(:,1),ones(1,size(WELL.Won,1)),QQmwo(:,1),na,1);
    QWOG.QQcwo(:,1) = QWOG.QQcwo(:,1) + sparse(WELL.WonC(:,1),ones(1,size(WELL.WonC,1)),QQcwo(:,1),nc,1);
    QWOG.QQgwo(:,1) = QWOG.QQgwo(:,1) + sparse(WELL.WonG(:,1),ones(1,size(WELL.WonG,1)),QQgwo(:,1),ng,1);
    QWOG.QQdwo(:,1) = QWOG.QQdwo(:,1) + sparse(WELL.WonD(:,1),ones(1,size(WELL.WonD,1)),QQdwo(:,1),nd,1);
   
    QWOG.QQmwo(:,2) = QWOG.QQmwo(:,2) + sparse(WELL.Won(:,1),ones(1,size(WELL.Won,1)),QQmwo(:,2),na,1);
    QWOG.QQcwo(:,2) = QWOG.QQcwo(:,2) + sparse(WELL.WonC(:,1),ones(1,size(WELL.WonC,1)),QQcwo(:,2),nc,1);
    QWOG.QQgwo(:,2) = QWOG.QQgwo(:,2) + sparse(WELL.WonG(:,1),ones(1,size(WELL.WonG,1)),QQgwo(:,2),ng,1);
    QWOG.QQdwo(:,2) = QWOG.QQdwo(:,2) + sparse(WELL.WonD(:,1),ones(1,size(WELL.WonD,1)),QQdwo(:,2),nd,1);
   
    Phi(1:na+nc+ng+nd,1) = Phi(1:na+nc+ng+nd,1) + dPt(1:na+nc+ng+nd);   %//////////////
    Phi(1:na+nc+ng+nd,2) = Phi(1:na+nc+ng+nd,2) + dPt(1:na+nc+ng+nd);   %//////////////
   %Phi(:,3) = Phi(:,3) + dPt(1:na+nc+ng+nd);
  
   [SGMMD]=SGim(GEOM.dV,Sw,So,SGM.Mp(:,2),PR.zc,SGM.Rs(:,2),PR.rs,SGM.Bwog,dPt,Pi,1,VEC.va,zeros(0,1),zeros(0,1),VEC.vd,VEC.vb,dt);
   Pi(1:na+nc+ng+nd,1) = Pi(1:na+nc+ng+nd,1) + dPt(1:na+nc+ng+nd);   
   if isempty(RC.Cr)~=0 || isempty(RC.Gr)~=0
     flag_gim = 0;
   else
     Fwater = AMW*Phi(1:na+nc+ng+nd,1) - dItimeW + [QWOG.QQmwo(:,1);QWOG.QQcwo(:,1);QWOG.QQgwo(:,1);QWOG.QQdwo(:,1)];
     Foil = AMO*Phi(1:na+nc+ng+nd,2) - dItimeO + [QWOG.QQmwo(:,2);QWOG.QQcwo(:,2);QWOG.QQgwo(:,2);QWOG.QQdwo(:,2)];
     Sw = Sw + (Fwater - SGM.Cwp.*dPt(1:na+nc+ng+nd))./SGM.Cws;
     So = 1 - Sw;
     Qm2 = QBild(QWOG.QQm(WELL.Won(:,3)),QWOG.QQmwo(WELL.Won(:,1),:),WELL.Uf(WELL.Won(:,3),ft+1),WELL.Won(:,1),dt,WELL.Won(:,3),nw,W1);
     Qd2 = QBild(QWOG.QQd(WELL.WonD(:,3)),QWOG.QQdwo(WELL.WonD(:,1),:),WELL.Uf(WELL.WonD(:,3),ft+1),WELL.WonD(:,1),dt,WELL.WonD(:,3),nw,W1D);  
     flag_gim=sum(abs(dPt(1:na+nc+ng+nd)./Pi(1:na+nc+ng+nd))>=1e-6)~=0;
   end;
  end; 
 
   
   
