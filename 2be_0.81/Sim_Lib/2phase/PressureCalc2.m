function [Pi,Sw,Pw,TL,DL,Phi,CMP,Qm2,Qc2,Qg2,Qd2,QQ,QQBND,QQoBND]=PressureCalc2(Pi,Sw,Phi,Sw0,Pw,Cp,TRM,KWOG,KWOG_GY,CMP,RC,WELL,fp,VEC,GEOM,DATA,M2FR,ft,PR,BXYZ,dt,Qf,GYData,BB,QQ,QQBND,QQoBND,ndt)
  na = RC.na;
  nc = RC.nc;
  ng = RC.ng;
  nd = RC.nd;
  nb = RC.nb;
  nw = size(Qf,1);
  Nsum = na+nc+ng+nd+nb;
  b1wb=sparse(nb,1);
  dPw = zeros(nw,1);
  dPt = zeros(Nsum,1);
  vad = RC.ADr;
  Qm2=zeros(nw,5);
  Qd2=zeros(nw,5);
  kj = 0;
  flag_gim = 1;
  while flag_gim==1 && kj<5
   kj=kj+1; 
   [SGM,CMP]=SGim2(GEOM.dV,Sw,PR.zc,Pi,dPt(1:Nsum),dt,CMP);
  
   [TW,TO,TP]=Potok_MKT2(TRM.TTM,Phi(VEC.va,:),KWOG,Cp(VEC.va,1),PR,RC.Arc2,fp,PR.kms(1),GEOM.L,CMP);  %проводимости по фазам
   [CW,CO,~]=Potok_Tube2(TRM.TC,Phi(VEC.vc,:),KWOG,Cp(VEC.vc),PR,fp,PR.kms(2),DATA.Lc,RC.Cr2,RC.Cc2,nc,CMP,VEC.vc);
   [GW,GO,~]=Potok_Tube2(TRM.TG,Phi(VEC.vg,:),KWOG,Cp(VEC.vg),PR,fp,PR.kms(3),GEOM.Lg,RC.Gr2,RC.Gc2,ng,CMP,VEC.vg);
   [DW,DO,DP]=Potok_Tube2(TRM.TD,Phi(VEC.vd,:),KWOG,Cp(VEC.vd),PR,fp,PR.kms(4),GEOM.Ld,RC.Dr2,RC.Dc2,nd,CMP,VEC.vd);
                                                   
   [A2CW,A2CO,~]=Obmen_T2M2(M2FR.A2C,Phi(VEC.va,:),Phi(VEC.vc,:),GEOM.K(:,1),PR,Cp(VEC.va),Cp(VEC.vc),KWOG,CMP,VEC.va,VEC.vc);
   [A2GW,A2GO,~]=Obmen_T2M2(M2FR.A2G,Phi(VEC.va,:),Phi(VEC.vg,:),GEOM.K(:,1),PR,Cp(VEC.va),Cp(VEC.vg),KWOG,CMP,VEC.va,VEC.vg);
        
   [A2DW,A2DO,A2DP]=Obmen_T2M2(M2FR.A2D,Phi(VEC.va,:),Phi(VEC.vd,:),GEOM.K(:,1),PR,Cp(VEC.va),Cp(VEC.vd),KWOG,CMP,VEC.va,VEC.vd);
   [A2BW,A2BO,A2BP]=Obmen_T2M2(M2FR.A2B,Phi(VEC.va,:),Phi(VEC.vb,:),ones(na,1),PR,Cp(VEC.va),Cp(VEC.vb),KWOG,CMP,VEC.va,VEC.vb);
   [D2BW,D2BO,D2BP]=Obmen_T2M2(M2FR.D2B,Phi(VEC.vd,:),Phi(VEC.vb,:),GEOM.K(:,1),PR,Cp(VEC.vd),Cp(VEC.vb),KWOG,CMP,VEC.vd,VEC.vb);
        
   [D2CW,D2CO,~]=Obmen_T2M2(M2FR.C2D,Phi(VEC.vd,:),Phi(VEC.vc,:),GEOM.K(:,1),PR,Cp(VEC.vd),Cp(VEC.vc),KWOG,CMP,VEC.vd,VEC.vc);
   [D2GW,D2GO,~]=Obmen_T2M2(M2FR.G2D,Phi(VEC.vd,:),Phi(VEC.vg,:),GEOM.K(:,1),PR,Cp(VEC.vd),Cp(VEC.vg),KWOG,CMP,VEC.vd,VEC.vg);
        
        %Wf=Wf0.*(1-0.0001*(P0(Won)-Pi(Won))).^3;
   [WBND,QQBND]=GY_bild2(GYData,Pi(:,1),dPt(:,1),KWOG_GY,Cp(:,1),RC,TRM.Txyz_GY_A,TRM.Txyz_GY_D,PR,CMP,VEC,QQBND);
      
   [W1,W6,Wo,W7]=Well_MKT2(WELL.Won,WELL.Uf(WELL.Won(:,3),ft),Cp(VEC.va),PR,WELL.CpW(WELL.Won(:,3),ft),CMP,KWOG,VEC.va);% W1 - проводимость по всей жидкости, W6 - только для воды, W7 - полимер 
   [W1C,W6C,WoC,W7C]=Well_MKT2(WELL.WonC,WELL.Uf(WELL.WonC(:,3),ft),Cp(VEC.vc),PR,WELL.CpW(WELL.WonC(:,3),ft),CMP,KWOG,VEC.vc);
   [W1G,W6G,WoG,W7G]=Well_MKT2(WELL.WonG,WELL.Uf(WELL.WonG(:,3),ft),Cp(VEC.vg),PR,WELL.CpW(WELL.WonG(:,3),ft),CMP,KWOG,VEC.vg);
   [W1D,W6D,WoD,W7D]=Well_MKT2(WELL.WonD,WELL.Uf(WELL.WonD(:,3),ft),Cp(VEC.vd),PR,WELL.CpW(WELL.WonD(:,3),ft),CMP,KWOG,VEC.vd);
        
   A1= -sparse(WELL.Won(:,1),WELL.Won(:,1),W1,na,na)-sparse(1:na,1:na,SGM.Clp(VEC.va)+sum(WBND.b1gm(:,1:2),2),na,na);  %Матрица коэф. для пор
   C1= -sparse(1:nc,1:nc,SGM.Clp(VEC.vc)',nc,nc)-sparse(WELL.WonC(:,1),WELL.WonC(:,1),W1C,nc,nc);                       %Матрица коэф. для вертикальных трещ.
   G1= -sparse(1:ng,1:ng,SGM.Clp(VEC.vg)',ng,ng)-sparse(WELL.WonG(:,1),WELL.WonG(:,1),W1G,ng,ng);                       %Матрица коэф. для гориз. трещ.
   D1= -sparse(WELL.WonD(:,1),WELL.WonD(:,1),W1D,nd,nd)-sparse(1:nd,1:nd,SGM.Clp(VEC.vd)+sum(WBND.b1gd(:,1:2),2),nd,nd);                       %Матрица коэф. для двойной пор.
   B1= -sparse(1:nb,1:nb,SGM.Clp(VEC.vb)'+WBND.b1gb',nb,nb);  %Матрица коэф. для границ
        
   AW1=TW-sparse(1:na,1:na,sum(A2CW,2)+sum(A2GW,2)+sum(A2BW,2)+sum(A2DW,2),na,na);  
   CW1=CW-sparse(1:nc,1:nc,sum(A2CW,1)+sum(D2CW,1),nc,nc);                       
   GW1=GW-sparse(1:ng,1:ng,sum(A2GW,1)+sum(D2GW,1),ng,ng);   
   DW1=DW-sparse(1:nd,1:nd,sum(A2DW,1)'+sum(D2BW,2)+sum(D2CW,2)+sum(D2GW,2),nd,nd);
   BW1=-sparse(1:nb,1:nb,sum(A2BW,1)+sum(D2BW,1),nb,nb);
   
   AO1=TO-sparse(1:na,1:na,sum(A2CO,2)+sum(A2GO,2)+sum(A2BO,2)+sum(A2DO,2),na,na);  
   CO1=CO-sparse(1:nc,1:nc,sum(A2CO,1)+sum(D2CO,1),nc,nc);                       
   GO1=GO-sparse(1:ng,1:ng,sum(A2GO,1)+sum(D2GO,1),ng,ng);   
   DO1=DO-sparse(1:nd,1:nd,sum(A2DO,1)'+sum(D2BO,2)+sum(D2CO,2)+sum(D2GO,2),nd,nd);
   
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
   
   Awater = sparse(1:Nsum,1:Nsum,CMP.Cw,Nsum,Nsum)*AMW;
            
   [dItime,dItimeW] = NR_Time2(GEOM.dV,Sw,Sw0,CMP,dt);         
   dPhi = Awater*Phi(:,1) + AMO*Phi(:,2);
 %  while  flag_pwq==1

    b1wm=sparse(WELL.Won(:,1),ones(1,size(WELL.Won,1)),-W1.*dPw(WELL.Won(:,3)),na,1);
    b1wc=sparse(WELL.WonC(:,1),ones(1,size(WELL.WonC,1)),-W1C.*dPw(WELL.WonC(:,3)),nc,1);
    b1wg=sparse(WELL.WonG(:,1),ones(1,size(WELL.WonG,1)),-W1G.*dPw(WELL.WonG(:,3)),ng,1);
    b1wd=sparse(WELL.WonD(:,1),ones(1,size(WELL.WonD,1)),-W1D.*dPw(WELL.WonD(:,3)),nd,1);
            
    b1wm=b1wm.*(sum(W2M(Qf==0,:),1)~=0)';
    b1wc=b1wc.*(sum(W2C(Qf==0,:),1)~=0)';
    b1wg=b1wg.*(sum(W2G(Qf==0,:),1)~=0)';
    b1wd=b1wd.*(sum(W2D(Qf==0,:),1)~=0)';
           
    QQ.QQm = QQ.QQm + sparse(WELL.Won(:,1),ones(1,size(W1,1)),W1.*(dPw(WELL.Won(:,3))-dPt(VEC.va(WELL.Won(:,1)))),na,1);
    QQ.QQc = QQ.QQc + sparse(WELL.WonC(:,1),ones(1,size(W1C,1)),W1C.*(dPw(WELL.WonC(:,3))-dPt(VEC.vc(WELL.WonC(:,1)))),nc,1);
    QQ.QQg = QQ.QQg + sparse(WELL.WonG(:,1),ones(1,size(W1G,1)),W1G.*(dPw(WELL.WonG(:,3))-dPt(VEC.vg(WELL.WonG(:,1)))),ng,1);
    QQ.QQd = QQ.QQd + sparse(WELL.WonD(:,1),ones(1,size(W1D,1)),W1D.*(dPw(WELL.WonD(:,3))-dPt(VEC.vd(WELL.WonD(:,1)))),nd,1);
   
    b1wm = b1wm - QQ.QQm;
    b1wc = b1wc - QQ.QQc;
    b1wg = b1wg - QQ.QQg;
    b1wd = b1wd - QQ.QQd;
                           
    BM = [b1wm;b1wc;b1wg;b1wd;b1wb] - [QQBND.Qm;zeros(nc,1);zeros(ng,1);QQBND.Qd;QQBND.Qb] - dPhi + dItime;
  
    W2M1=WM1(Qf~=0,:);
    W2M2=WM2(:,Qf~=0);
    W2M3=WM3(Qf~=0,Qf~=0);
                  
    AM = [A1, zeros(na,nc), zeros(na,ng), zeros(na,nd), zeros(na,nb);
          zeros(nc,na), C1, zeros(nc,ng), zeros(nc,nd), zeros(nc,nb);
          zeros(ng,na), zeros(ng,nc), G1, zeros(ng,nd), zeros(ng,nb);
          zeros(nd,na), zeros(nd,nc), zeros(nd,ng), D1, zeros(nd,nb);
          zeros(nb,na), zeros(nb,nc), zeros(nb,ng), zeros(nb,nd), B1+BB];
      
    AM = AM + Awater + AMO; 
    
    Qz1 = sparse(WELL.Won(Qf~=0,3),ones(1,size(WELL.Won(Qf~=0),1)),Qf(Qf~=0),nw,1);
    Qm1 = sparse(WELL.Won(Qf~=0,3),ones(1,size(WELL.Won(Qf~=0),1)),QQ.QQm(WELL.Won(Qf~=0),1),nw,1);
    Qc1 = sparse(WELL.WonC(Qf(WELL.WonC(:,3))~=0,3),ones(1,size(WELL.WonC(Qf(WELL.WonC(:,3))~=0,3),1)),QQ.QQc(WELL.WonC(Qf(WELL.WonC(:,3))~=0,1)),nw,1);
    Qg1 = sparse(WELL.WonG(Qf(WELL.WonG(:,3))~=0,3),ones(1,size(WELL.WonG(Qf(WELL.WonG(:,3))~=0,3),1)),QQ.QQg(WELL.WonG(Qf(WELL.WonG(:,3))~=0,1)),nw,1);
    Qd1 = sparse(WELL.WonD(Qf(WELL.WonD(:,3))~=0,3),ones(1,size(WELL.WonD(Qf(WELL.WonD(:,3))~=0,3),1)),QQ.QQd(WELL.WonD(Qf(WELL.WonD(:,3))~=0,1)),nw,1);
   
    Q1 = -Qz1+Qm1+Qc1+Qg1+Qd1;
    dPt = [BM',Q1(Qf~=0)']/[AM,W2M2;W2M1,W2M3];  %/////////////
    dPt = dPt';
      
    dPw(Qf~=0) = dPt(Nsum+1:end);
   
    [QQoBND] = QIterGY2(dPt,QQoBND,WBND,BXYZ);
    
    [QQ.QQmwo] = QIter2(QQ.QQmwo,W6,Wo,dPt(VEC.va),WELL.Won(:,1),dPw(WELL.Won(:,3)),na); 
    [QQ.QQcwo] = QIter2(QQ.QQcwo,W6C,WoC,dPt(VEC.vc),WELL.WonC(:,1),dPw(WELL.WonC(:,3)),nc);
    [QQ.QQgwo] = QIter2(QQ.QQgwo,W6G,WoG,dPt(VEC.vg),WELL.WonG(:,1),dPw(WELL.WonG(:,3)),ng);
    [QQ.QQdwo] = QIter2(QQ.QQdwo,W6D,WoD,dPt(VEC.vd),WELL.WonD(:,1),dPw(WELL.WonD(:,3)),nd);
       
    Phi(1:Nsum,1) = Phi(1:Nsum,1) + dPt(1:Nsum); 
    Phi(1:Nsum,2) = Phi(1:Nsum,2) + dPt(1:Nsum); 
    Pi(1:Nsum,1) = Pi(1:Nsum,1) + dPt(1:Nsum);  
    Pw = Pw + dPw;  
       
    Fwater = AMW*Phi(1:Nsum,1) - dItimeW + [QQ.QQmwo(:,1);QQ.QQcwo(:,1);QQ.QQgwo(:,1);QQ.QQdwo(:,1);zeros(nb,1)] + [QQoBND.Qmw;zeros(nc,1);zeros(ng,1);QQoBND.Qdw;zeros(nb,1)];
   
    Sw = Sw + (Fwater - SGM.Cwp.*dPt(1:Nsum))./SGM.Cwsw;
         
    flag_gim=sum(abs(dPt(1:na+nc+ng+nd)./Pi(1:na+nc+ng+nd))>=1e-6)~=0; 
   end;
  %end; 
   Qm2 = QBild(QQ.QQm(WELL.Won(:,1)),QQ.QQmwo(WELL.Won(:,1),:),WELL.Uf(WELL.Won(:,3),ft),WELL.Won(:,1),dt,WELL.Won(:,3),nw,W1);
   Qc2 = QBild(QQ.QQc(WELL.WonC(:,1)),QQ.QQcwo(WELL.WonC(:,1),:),WELL.Uf(WELL.WonC(:,3),ft),WELL.WonC(:,1),dt,WELL.WonC(:,3),nw,W1C);
   Qg2 = QBild(QQ.QQg(WELL.WonG(:,1)),QQ.QQgwo(WELL.WonG(:,1),:),WELL.Uf(WELL.WonG(:,3),ft),WELL.WonG(:,1),dt,WELL.WonG(:,3),nw,W1G);
   Qd2 = QBild(QQ.QQd(WELL.WonD(:,1)),QQ.QQdwo(WELL.WonD(:,1),:),WELL.Uf(WELL.WonD(:,3),ft),WELL.WonD(:,1),dt,WELL.WonD(:,3),nw,W1D);
  TL = sparse(1:na,1:na,CMP.Cw(VEC.va),na,na)*TW + TO;
  DL = sparse(1:nd,1:nd,CMP.Cw(VEC.vd),nd,nd)*DW + DO;
  
     
 
   
   
