function [Pi,Pb,Sw,So,CL,GL,Phi,SGM,CMP,Qm2,Qd2,VSAT,QQ,QQBND,QQwoBND]=PressureCalc2(Pi,Sw,So,Phi,Sw0,So0,Pb,Pw,Cp,TRM,KWOG,KWOG_GY,SGM,CMP,RC,WELL,fp,VEC,GEOM,DATA,NCELL,M2FR,ft,PR,BXYZ,dt,t,Qf,Qz,VSAT,GYData,BB,QQ,QQBND,QQwoBND)

  na = NCELL.na;
  nc = NCELL.nc;
  ng = NCELL.ng;
  nd = NCELL.nd;
  nb = NCELL.nb;
  nw = NCELL.nw;
  Nsum = na+nc+ng+nd+nb;
  Uf = WELL.Uf(:,ft+1);
  CpW = WELL.CpW(:,ft+1);
  b1wb=sparse(nb,1);
  vad = RC.ADr;
  dPw = zeros(size(Qz(:,1),1),1);
  dPt = zeros(Nsum,1);
  Qm2=zeros(nw,5);
  Qd2=zeros(nw,5);
  kj = 0;
  flag_gim = 1;
  while flag_gim==1 && kj<20
   kj=kj+1;
   Pbl = Pb;
   Sol = So;
   Swl = Sw;
   
   [TW,TO,TGG,TORS,TP]=Potok_MKT(TRM.TTM,Phi(VEC.va,:),KWOG,Cp(VEC.va,1),PR,RC.Arc2,fp,PR.kms(1),GEOM.L,SGM,CMP);  %проводимости по фазам
   [CW,CO,CG,CWCgs,COCgs,~]=Potok_Tube(TRM.TC,Phi(VEC.vc,:),Sw(VEC.vc),Cp(VEC.vc),PR,fp,PR.kms(2),DATA.Lc,RC.Cr2,RC.Cc2,NCELL.nc,SGM,CMP);
   [GW,GO,GG,GWCgs,GOCgs,~]=Potok_Tube(TRM.TG,Phi(VEC.vg,:),Sw(VEC.vg),Cp(VEC.vg),PR,fp,PR.kms(3),GEOM.Lg,RC.Gr2,RC.Gc2,NCELL.ng,SGM,CMP);
   [DW,DO,DG,DWCgs,DOCgs,DP]=Potok_Tube(TRM.TD,Phi(VEC.vd,:),Sw(VEC.vd),Cp(VEC.vd),PR,fp,PR.kms(4),GEOM.Ld,RC.Dr2,RC.Dc2,NCELL.nd,SGM,CMP);
                                                   
   [A2CW,A2CO,A2CG,A2CWCgs,A2COCgs,~]=Obmen_T2M(M2FR.A2C,Pi(VEC.va),Pi(VEC.vc),Sw(VEC.va),Sw(VEC.vc),GEOM.K(:,1),PR,Cp(VEC.va),Cp(VEC.vc),SGM,CMP);
   [A2GW,A2GO,A2GG,A2GWCgs,A2GOCgs,~]=Obmen_T2M(M2FR.A2G,Pi(VEC.va),Pi(VEC.vg),Sw(VEC.va),Sw(VEC.vg),GEOM.K(:,1),PR,Cp(VEC.va),Cp(VEC.vg),SGM,CMP);
        
   [A2DW,A2DO,A2DG,A2DWCgs,A2DOCgs,A2DP]=Obmen_T2M(M2FR.A2D,Pi(VEC.va),Pi(VEC.vd),Sw(VEC.va),Sw(VEC.vd),GEOM.K(:,1),PR,Cp(VEC.va),Cp(VEC.vd),SGM,CMP);
   [A2BW,A2BO,A2BG,A2BWCgs,A2BOCgs,A2BP]=Obmen_T2M(M2FR.A2B,Pi(VEC.va),Pi(VEC.vb),Sw(VEC.va),Sw(VEC.vb),ones(na,1),PR,Cp(VEC.va),Cp(VEC.vb),SGM,CMP);
   [D2BW,D2BO,D2BG,D2BWCgs,D2BOCgs,D2BP]=Obmen_T2M(M2FR.D2B,Pi(VEC.vd),Pi(VEC.vb),Sw(VEC.vd),Sw(VEC.vb),GEOM.K(:,1),PR,Cp(VEC.vd),Cp(VEC.vb),SGM,CMP);
        
   [D2CW,D2CO,D2CG,D2CWCgs,D2COCgs,~]=Obmen_T2M(M2FR.C2D,Pi(VEC.vd),Pi(VEC.vc),Sw(VEC.vd),Sw(VEC.vc),GEOM.K(:,1),PR,Cp(VEC.vd),Cp(VEC.vc),SGM,CMP);
   [D2GW,D2GO,D2GG,D2GWCgs,D2GOCgs,~]=Obmen_T2M(M2FR.G2D,Pi(VEC.vd),Pi(VEC.vg),Sw(VEC.vd),Sw(VEC.vg),GEOM.K(:,1),PR,Cp(VEC.vd),Cp(VEC.vg),SGM,CMP);
        
        %Wf=Wf0.*(1-0.0001*(P0(Won)-Pi(Won))).^3;
   [WBND,QQBND]=GY_bild(GYData,Pi([VEC.va,VEC.vd],1),dPt([VEC.va,VEC.vd],1),KWOG_GY,Cp([VEC.va,VEC.vd],1),RC,TRM.Txyz_GY_A,TRM.Txyz_GY_D,PR,SGM,CMP,VEC,QQBND);
      
   [W1,W6,Wo,Wg,W7,QQ.QQm]=Well_MKT(dPt(VEC.va),dPw(WELL.Won(:,3)),WELL.Won,Uf(WELL.Won(:,3)),Sw(VEC.va),So(VEC.va),Cp(VEC.va),PR.aw,PR.as,PR,CpW(WELL.Won(:,3)),SGM,CMP,KWOG,QQ.QQm);% W1 - проводимость по всей жидкости, W6 - только для воды, W7 - полимер 
   [W1C,W6C,WoC,WgC,W7C,QQ.QQc]=Well_MKT(dPt(VEC.vc),dPw(WELL.WonC(:,3)),WELL.WonC,Uf(WELL.WonC(:,3)),Sw(VEC.vc),So(VEC.vc),Cp(VEC.vc),PR.tw,PR.ts,PR,CpW(WELL.WonC(:,3)),SGM,CMP,KWOG,QQ.QQc);
   [W1G,W6G,WoG,WgG,W7G,QQ.QQg]=Well_MKT(dPt(VEC.vg),dPw(WELL.WonG(:,3)),WELL.WonG,Uf(WELL.WonG(:,3)),Sw(VEC.vg),So(VEC.vg),Cp(VEC.vg),PR.tw,PR.ts,PR,CpW(WELL.WonG(:,3)),SGM,CMP,KWOG,QQ.QQg);
   [W1D,W6D,WoD,WgD,W7D,QQ.QQd]=Well_MKT(dPt(VEC.vd),dPw(WELL.WonD(:,3)),WELL.WonD,Uf(WELL.WonD(:,3)),Sw(VEC.vd),So(VEC.vd),Cp(VEC.vd),PR.tw,PR.ts,PR,CpW(WELL.WonD(:,3)),SGM,CMP,KWOG,QQ.QQd);
        
   A1= -sparse(WELL.Won(:,1),WELL.Won(:,1),W1 - PR.rs*SGM.Cg(WELL.Won(:,1)).*QQ.QQmwo(:,2),na,na)-sparse(1:na,1:na,SGM.Clp(VEC.va)+sum(WBND.b1gm(:,1:2),2),na,na);  %Матрица коэф. для пор
   C1= -sparse(1:nc,1:nc,SGM.Clp(VEC.vc)',nc,nc)-sparse(WELL.WonC(:,1),WELL.WonC(:,1),W1C,nc,nc);                       %Матрица коэф. для вертикальных трещ.
   G1= -sparse(1:ng,1:ng,SGM.Clp(VEC.vg)',ng,ng)-sparse(WELL.WonG(:,1),WELL.WonG(:,1),W1G,ng,ng);                       %Матрица коэф. для гориз. трещ.
   D1= -sparse(WELL.WonD(:,1),WELL.WonD(:,1),W1D,nd,nd)-sparse(1:nd,1:nd,SGM.Clp(VEC.vd)+sum(WBND.b1gd(:,1:2),2),nd,nd);                       %Матрица коэф. для двойной пор.
   B1= -sparse(1:nb,1:nb,SGM.Clp(VEC.vb)'+WBND.b1gb',nb,nb);  %Матрица коэф. для границ
        
   AG1=TGG-sparse(1:na,1:na,sum(A2CG,2)+sum(A2GG,2)+sum(A2BG,2)+sum(A2DG,2),na,na);  
   CG1=CG-sparse(1:nc,1:nc,sum(A2CG,1)+sum(D2CG,1),nc,nc);                       
   GG1=GG-sparse(1:ng,1:ng,sum(A2GG,1)+sum(D2GG,1),ng,ng);   
   DG1=DG-sparse(1:nd,1:nd,sum(A2DG,1)'+sum(D2BG,2)+sum(D2CG,2)+sum(D2GG,2),nd,nd);
   
   AW1=TW-sparse(1:na,1:na,sum(A2CW,2)+sum(A2GW,2)+sum(A2BW,2)+sum(A2DW,2),na,na);  
   CW1=CW-sparse(1:nc,1:nc,sum(A2CW,1)+sum(D2CW,1),nc,nc);                       
   GW1=GW-sparse(1:ng,1:ng,sum(A2GW,1)+sum(D2GW,1),ng,ng);   
   DW1=DW-sparse(1:nd,1:nd,sum(A2DW,1)'+sum(D2BW,2)+sum(D2CW,2)+sum(D2GW,2),nd,nd);
   BW1=-sparse(1:nb,1:nb,sum(A2BW,1)+sum(D2BW,1),nb,nb);
   
   AO1=TO-sparse(1:na,1:na,sum(A2CO,2)+sum(A2GO,2)+sum(A2BO,2)+sum(A2DO,2),na,na);  
   CO1=CO-sparse(1:nc,1:nc,sum(A2CO,1)+sum(D2CO,1),nc,nc);                       
   GO1=GO-sparse(1:ng,1:ng,sum(A2GO,1)+sum(D2GO,1),ng,ng);   
   DO1=DO-sparse(1:nd,1:nd,sum(A2DO,1)'+sum(D2BO,2)+sum(D2CO,2)+sum(D2GO,2),nd,nd);
   
   AMG=[AG1,   A2CG, A2GG, A2DG, A2BG;
       A2CG',  CG1,  M2FR.C2G, D2CG', M2FR.C2B;
       A2GG', M2FR.C2G', GG1,  D2GG', M2FR.G2B;
       A2DG', D2CG, D2GG,  DG1,  D2BG;
       A2BG', M2FR.C2B',M2FR.G2B',D2BG', zeros(nb,nb)];
   
   AMW=[AW1,   A2CW, A2GW, A2DW, A2BW;
       A2CW',  CW1,  M2FR.C2G, D2CW', M2FR.C2B;
       A2GW', M2FR.C2G', GW1,  D2GW', M2FR.G2B;
       A2DW', D2CW, D2GW,  DW1,  D2BW;
       A2BW', M2FR.C2B',M2FR.G2B',D2BW', BW1];
   
   AMO=[AO1,   A2CO, A2GO, A2DO, A2BO;
       A2CO',  CO1,  M2FR.C2G, D2CO', M2FR.C2B;
       A2GO', M2FR.C2G', GO1,  D2GO', M2FR.G2B;
       A2DO', D2CO, D2GO,  DO1,  D2BO;
       A2BO', M2FR.C2B',M2FR.G2B',D2BO', zeros(nb,nb)];
   
   AMORS=[TORS,   A2CO, A2GO, A2DO, A2BO;
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
   
    Aw = sparse(1:Nsum,1:Nsum,SGM.Cw,Nsum,Nsum)*AMW;
    Aoil = sparse(1:Nsum,1:Nsum,SGM.Co,Nsum,Nsum)*AMO;
    Aoilrs = sparse(1:Nsum,1:Nsum,SGM.Cor,Nsum,Nsum)*AMORS;
    Agas = sparse(1:Nsum,1:Nsum,SGM.Cg,Nsum,Nsum)*AMG;        
            
   [dItime,dItimeW,dItimeO,dItimeORs] = NR_Time(GEOM.dV,Sw,So,Sw0,So0,SGM,CMP,dt);         
   dPhi = Aw*Phi(:,1) + (Aoil + Aoilrs)*Phi(:,2) + Agas*Phi(:,3);
 %  while  flag_pwq==1

    b1wm=sparse(WELL.Won(:,1),ones(1,size(WELL.Won,1)),-W1.*dPw(WELL.Won(:,3)),na,1);
    b1wc=sparse(WELL.WonC(:,1),ones(1,size(WELL.WonC,1)),-W1C.*dPw(WELL.WonC(:,3)),nc,1);
    b1wg=sparse(WELL.WonG(:,1),ones(1,size(WELL.WonG,1)),-W1G.*dPw(WELL.WonG(:,3)),ng,1);
    b1wd=sparse(WELL.WonD(:,1),ones(1,size(WELL.WonD,1)),-W1D.*dPw(WELL.WonD(:,3)),nd,1);
            
    b1wm=b1wm.*(sum(W2M(Qf==0,:),1)~=0)';
    b1wc=b1wc.*(sum(W2C(Qf==0,:),1)~=0)';
    b1wg=b1wg.*(sum(W2G(Qf==0,:),1)~=0)';
    b1wd=b1wd.*(sum(W2D(Qf==0,:),1)~=0)';
            
    b1wm = b1wm - sparse(WELL.Won(:,1),ones(1,size(WELL.Won,1)),QQ.QQm,na,1);
    b1wc = b1wc - sparse(WELL.WonC(:,1),ones(1,size(WELL.WonC,1)),QQ.QQc,nc,1);
    b1wg = b1wg - sparse(WELL.WonG(:,1),ones(1,size(WELL.WonG,1)),QQ.QQg,ng,1);
    b1wd = b1wd - sparse(WELL.WonD(:,1),ones(1,size(WELL.WonD,1)),QQ.QQd,nd,1);
                           
    BM = [b1wm;b1wc;b1wg;b1wd;b1wb] - [QQBND.Qm;zeros(nc,1);zeros(ng,1);QQBND.Qd;QQBND.Qb] - dPhi + dItime;
  
    W2M1=WM1(Qf~=0,:);
    W2M2=WM2(:,Qf~=0);
    W2M3=WM3(Qf~=0,Qf~=0);
                  
    AM = [A1, zeros(na,nc), zeros(na,ng), zeros(na,nd), zeros(na,nb);
          zeros(nc,na), C1, zeros(nc,ng), zeros(nc,nd), zeros(nc,nb);
          zeros(ng,na), zeros(ng,nc), G1, zeros(ng,nd), zeros(ng,nb);
          zeros(nd,na), zeros(nd,nc), zeros(nd,ng), D1, zeros(nd,nb);
          zeros(nb,na), zeros(nb,nc), zeros(nb,ng), zeros(nb,nd), B1+BB];
      
    AM = AM + Aw + Aoil + Aoilrs + Agas; 
    
    Qz1 = sparse(WELL.Won(Qf~=0,3),ones(1,size(WELL.Won(Qf~=0),1)),Qz(Qf~=0,ft+1),nw,1);
    Qm1 = sparse(WELL.Won(Qf~=0,3),ones(1,size(WELL.Won(Qf~=0),1)),QQ.QQm(Qf~=0),nw,1);
    Qc1 = sparse(WELL.WonC(Qf(WELL.WonC(:,3))~=0,3),ones(1,size(WELL.WonC(Qf(WELL.WonC(:,3))~=0,3),1)),QQ.QQc(Qf(WELL.WonC(:,3))~=0),nw,1);
    Qg1 = sparse(WELL.WonG(Qf(WELL.WonG(:,3))~=0,3),ones(1,size(WELL.WonG(Qf(WELL.WonG(:,3))~=0,3),1)),QQ.QQg(Qf(WELL.WonG(:,3))~=0),nw,1);
    Qd1 = sparse(WELL.WonD(Qf(WELL.WonD(:,3))~=0,3),ones(1,size(WELL.WonD(Qf(WELL.WonD(:,3))~=0,3),1)),QQ.QQd(Qf(WELL.WonD(:,3))~=0),nw,1);
   
    Q1 = -Qz1+Qm1+Qc1+Qg1+Qd1;
    dPt = [BM',Q1(Qf~=0)']/[AM,W2M2;W2M1,W2M3];  %/////////////
    dPt = dPt';
      
    dPw(Qf~=0) = dPt(Nsum+1:end);
   
    [QQwoBND] = QIterGY(dPt,QQwoBND,WBND,SGM,CMP,BXYZ);
    
    [QQ.QQmwo] = QIter(QQ.QQmwo,W6,Wo,dPt(VEC.va),WELL.Won(:,1),dPw(WELL.Won(:,3)),SGM,CMP); 

    QQmwog(:,1) =  sparse(WELL.Won(:,1),ones(1,size(WELL.Won,1)),QQ.QQmwo(:,1),na,1);
    QQcwog(:,1) =  sparse(WELL.WonC(:,1),ones(1,size(WELL.WonC,1)),QQ.QQcwo(:,1),nc,1);
    QQgwog(:,1) =  sparse(WELL.WonG(:,1),ones(1,size(WELL.WonG,1)),QQ.QQgwo(:,1),ng,1);
    QQdwog(:,1) =  sparse(WELL.WonD(:,1),ones(1,size(WELL.WonD,1)),QQ.QQdwo(:,1),nd,1);
    
    QQmwog(:,2) =  sparse(WELL.Won(:,1),ones(1,size(WELL.Won,1)),QQ.QQmwo(:,2),na,1);
    QQcwog(:,2) =  sparse(WELL.WonC(:,1),ones(1,size(WELL.WonC,1)),QQ.QQcwo(:,2),nc,1);
    QQgwog(:,2) =  sparse(WELL.WonG(:,1),ones(1,size(WELL.WonG,1)),QQ.QQgwo(:,2),ng,1);
    QQdwog(:,2) =  sparse(WELL.WonD(:,1),ones(1,size(WELL.WonD,1)),QQ.QQdwo(:,2),nd,1);  
   
 %  if isempty(RC.Cr)~=0 || isempty(RC.Gr)~=0
 %    flag_gim = 0;
 %  else
    Phi(1:Nsum,1) = Phi(1:Nsum,1) + dPt(1:Nsum); 
    Phi(1:Nsum,2) = Phi(1:Nsum,2) + dPt(1:Nsum); 
    Phi(1:Nsum,3) = Phi(1:Nsum,3) + dPt(1:Nsum);
    Pi(1:Nsum,1) = Pi(1:Nsum,1) + dPt(1:Nsum); 
    
    Fwater = AMW*Phi(1:Nsum,1) - dItimeW + [QQmwog(:,1);QQcwog(:,1);QQgwog(:,1);QQdwog(:,1);zeros(nb,1)] + [QQwoBND.Qmw;zeros(nc,1);zeros(ng,1);QQwoBND.Qdw;QQBND.Qb - WBND.b1gb.*dPt(VEC.vb)];
    FoilRs = AMORS*Phi(1:Nsum,2) - dItimeORs + CMP.Rs(:,2).*([QQmwog(:,2);QQcwog(:,2);QQgwog(:,2);QQdwog(:,2);zeros(nb,1)] + [QQwoBND.Qmo;zeros(nc,1);zeros(ng,1);QQwoBND.Qdo;QQBND.Qb - WBND.b1gb.*dPt(VEC.vb)]);
    Foil = AMO*Phi(1:Nsum,2) - dItimeO + [QQmwog(:,2);QQcwog(:,2);QQgwog(:,2);QQdwog(:,2);zeros(nb,1)] + [QQwoBND.Qmo;zeros(nc,1);zeros(ng,1);QQwoBND.Qdo;zeros(nb,1)];
     
     dSw(VSAT.vp,1) = (Fwater(VSAT.vp) - SGM.Cwp(VSAT.vp).*dPt(VSAT.vp))./SGM.Cwsw(VSAT.vp);
     dSo(VSAT.vp,1) = (Foil(VSAT.vp) - SGM.Cop(VSAT.vp).*dPt(VSAT.vp))./SGM.Coso(VSAT.vp);
     dPb(VSAT.vp,1) = 0;
     
     dSw(VSAT.vg,1) = (Foil(VSAT.vg) - SGM.Cop(VSAT.vg).*dPt(VSAT.vg)).*SGM.Cgpb(VSAT.vg) - (FoilRs(VSAT.vg) - SGM.Cgp(VSAT.vg).*dPt(VSAT.vg)).*SGM.Copb(VSAT.vg);
     dPb(VSAT.vg,1) = (FoilRs(VSAT.vg) - SGM.Cgp(VSAT.vg).*dPt(VSAT.vg)).*SGM.Cosw(VSAT.vg) - (Foil(VSAT.vg) - SGM.Cop(VSAT.vg).*dPt(VSAT.vg)).*SGM.Cgsw(VSAT.vg);
     dSo(VSAT.vg,1) = 0;
     
     vg1 = Pbl + dPb > Pi + 1.e-10;
     vg = Pbl(VSAT.vg) + dPb(VSAT.vg) > Pi(VSAT.vg) + 1.e-10;     
     vp1 = Sol + Swl > 1 + 1.e-1;   
  %   vp = So(VSAT.vp) + dSo(VSAT.vp) + Sw(VSAT.vp) + dSw(VSAT.vp) > 1 + 1.e-3;
     vp = Sol(VSAT.vp) + Swl(VSAT.vp) > 1 + 1.e-1;  
     VSAT.vg(vg) = [];
     VSAT.vp(vp) = []; 
     VSAT.vp = [VSAT.vp,VEC.v(vg1)]; VSAT.vg = [VSAT.vg,VEC.v(vp1)];
     
     dPb(VEC.v(vp1)) = - 0.0001;
     dSo(VEC.v(vg1)) = - 0.0001;

     Sw(1:Nsum,1) = Swl(1:Nsum,1) + dSw(1:Nsum);   
     So(VSAT.vg,1) = 1 - dSw(VSAT.vg);    
     Pb(VSAT.vg) = Pbl(VSAT.vg) + dPb(VSAT.vg);
       
     So(VSAT.vp,1) = Sol(VSAT.vp) + dSo(VSAT.vp);    
     Pb(VSAT.vp) = Pi(VSAT.vp);
     
     Qm2 = QBild(QQ.QQm(WELL.Won(:,3)),QQmwog(WELL.Won(:,1),:),WELL.Uf(WELL.Won(:,3),ft+1),WELL.Won(:,1),dt,WELL.Won(:,3),nw,W1);
     Qd2 = QBild(QQ.QQd(WELL.WonD(:,3)),QQdwog(WELL.WonD(:,1),:),WELL.Uf(WELL.WonD(:,3),ft+1),WELL.WonD(:,1),dt,WELL.WonD(:,3),nw,W1D);  
     flag_gim=sum(abs(dPt(1:na+nc+ng+nd)./Pi(1:na+nc+ng+nd))>=1e-6)~=0;
     [SGM,CMP]=SGim(GEOM.dV,Sw,So,PR.zc,PR.rs,Pi,Pb,dPb,dPt,dt,VSAT,CMP);
  % end;
  end; 
  CL = CG + CWCgs + COCgs;
  GL = GG + GWCgs + GOCgs;
  
     
 
   
   
