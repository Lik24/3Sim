function [Pi,Pb,Sw,So,Pw,TL,DL,Phi,CMP,Qm2,Qc2,Qg2,Qd2,VSAT,QQ,QQBND,QQwoBND,kj]=PressureCalcBO(Pi,Sw,So,Phi,Sw0,So0,Pb,Pw,Cp,TRM,KWOG,KWOG_GY,CMP,RC,WELL,fp,VEC,GEOM,DATA,NCELL,M2FR,ft,PR,BXYZ,dt,Qf,VSAT,GYData,BB,QQ,QQBND,QQwoBND)
  na = NCELL.na;
  nc = NCELL.nc;
  ng = NCELL.ng;
  nd = NCELL.nd;
  nb = NCELL.nb;
  nw = NCELL.nw;
  Nsum = na+nc+ng+nd+nb;
  b1wb=sparse(nb,1);
  vad = RC.ADr;
  dPw = zeros(size(Qf,1),1);
  dPt = zeros(Nsum,1);
  dPb = zeros(Nsum,1);
  Qm2=zeros(nw,5);
  Qd2=zeros(nw,5);
  kj = 0;
  flag_gim = 1;
  while flag_gim==1 && kj<50
   kj=kj+1;
   Pbl = Pb;
   Sol = So;
   Swl = Sw;
   [SGM,CMP]=SGimBO(GEOM.dV,Sw,So,PR.zc,PR.rs,Pi,Pb,Pb-Pbl,dPt(1:Nsum),dt,VSAT,CMP);
   
   [TW,TO,TGG,TORS,TP]=Potok_MKTBO(TRM.TTM,Phi(VEC.va,:),KWOG,Cp(VEC.va,1),PR,RC.Arc2,fp,PR.kms(1),GEOM.L,CMP);  %проводимости по фазам
   [CW,CO,CG,CORS,~]=Potok_TubeBO(TRM.TC,Phi(VEC.vc,:),KWOG,Cp(VEC.vc),PR,fp,PR.kms(2),DATA.Lc,RC.Cr2,RC.Cc2,nc,CMP,VEC.vc);
   [GW,GO,GG,GORS,~]=Potok_TubeBO(TRM.TG,Phi(VEC.vg,:),KWOG,Cp(VEC.vg),PR,fp,PR.kms(3),GEOM.Lg,RC.Gr2,RC.Gc2,ng,CMP,VEC.vg);
   [DW,DO,DG,DORS,DP]=Potok_TubeBO(TRM.TD,Phi(VEC.vd,:),KWOG,Cp(VEC.vd),PR,fp,PR.kms(4),GEOM.Ld,RC.Dr2,RC.Dc2,nd,CMP,VEC.vd);
                                                   
   [A2CW,A2CO,A2CG,A2CORS,~]=Obmen_T2MBO(M2FR.A2C,Phi(VEC.va,:),Phi(VEC.vc,:),GEOM.K(:,1),PR,Cp(VEC.va),Cp(VEC.vc),KWOG,CMP,VEC.va,VEC.vc);
   [A2GW,A2GO,A2GG,A2GORS,~]=Obmen_T2MBO(M2FR.A2G,Phi(VEC.va,:),Phi(VEC.vg,:),GEOM.K(:,1),PR,Cp(VEC.va),Cp(VEC.vg),KWOG,CMP,VEC.va,VEC.vg);
      
   [A2DW,A2DO,A2DG,A2DORS,A2DP]=Obmen_T2MBO(M2FR.A2D,Phi(VEC.va,:),Phi(VEC.vd,:),GEOM.K(:,1),PR,Cp(VEC.va),Cp(VEC.vd),KWOG,CMP,VEC.va,VEC.vd);
   [A2BW,A2BO,A2BG,A2BORS,A2BP]=Obmen_T2MBO(M2FR.A2B,Phi(VEC.va,:),Phi(VEC.vb,:),ones(na,1),PR,Cp(VEC.va),Cp(VEC.vb),KWOG,CMP,VEC.va,VEC.vb);
   [D2BW,D2BO,D2BG,D2BORS,D2BP]=Obmen_T2MBO(M2FR.D2B,Phi(VEC.vd,:),Phi(VEC.vb,:),GEOM.K(:,1),PR,Cp(VEC.vd),Cp(VEC.vb),KWOG,CMP,VEC.vd,VEC.vb);
        
   [D2CW,D2CO,D2CG,D2CORS,~]=Obmen_T2MBO(M2FR.C2D,Phi(VEC.vd,:),Phi(VEC.vc,:),GEOM.K(:,1),PR,Cp(VEC.vd),Cp(VEC.vc),KWOG,CMP,VEC.vd,VEC.vc);
   [D2GW,D2GO,D2GG,D2GORS,~]=Obmen_T2MBO(M2FR.G2D,Phi(VEC.vd,:),Phi(VEC.vg,:),GEOM.K(:,1),PR,Cp(VEC.vd),Cp(VEC.vg),KWOG,CMP,VEC.vd,VEC.vg);
        
        %Wf=Wf0.*(1-0.0001*(P0(Won)-Pi(Won))).^3;
   [WBND,QQBND]=GY_bildBO(GYData,Pi([VEC.va,VEC.vd],1),dPt([VEC.va,VEC.vd],1),KWOG_GY,Cp([VEC.va,VEC.vd],1),RC,TRM.Txyz_GY_A,TRM.Txyz_GY_D,PR,CMP,VEC,QQBND);
      
   [W1,W6,Wo,Wg,W7]=Well_MKTBO(WELL.Won,WELL.Uf(WELL.Won(:,3),ft),Cp(VEC.va),PR,WELL.CpW(WELL.Won(:,3),ft),CMP,KWOG,VEC.va);% W1 - проводимость по всей жидкости, W6 - только для воды, W7 - полимер 
   [W1C,W6C,WoC,WgC,W7C]=Well_MKTBO(WELL.WonC,WELL.Uf(WELL.WonC(:,3),ft),Cp(VEC.vc),PR,WELL.CpW(WELL.WonC(:,3),ft),CMP,KWOG,VEC.vc);
   [W1G,W6G,WoG,WgG,W7G]=Well_MKTBO(WELL.WonG,WELL.Uf(WELL.WonG(:,3),ft),Cp(VEC.vg),PR,WELL.CpW(WELL.WonG(:,3),ft),CMP,KWOG,VEC.vg);
   [W1D,W6D,WoD,WgD,W7D]=Well_MKTBO(WELL.WonD,WELL.Uf(WELL.WonD(:,3),ft),Cp(VEC.vd),PR,WELL.CpW(WELL.WonD(:,3),ft),CMP,KWOG,VEC.vd);
        
   A1= -sparse(WELL.Won(:,1),WELL.Won(:,1),W1,na,na)-sparse(1:na,1:na,SGM.Clp(VEC.va)+sum(WBND.b1gm(:,1:2),2),na,na);
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
   
   AORS1=TORS-sparse(1:na,1:na,sum(A2CORS,2)+sum(A2GORS,2)+sum(A2BORS,2)+sum(A2DORS,2),na,na);  
   CORS1=CORS-sparse(1:nc,1:nc,sum(A2CORS,1)+sum(D2CORS,1),nc,nc);                       
   GORS1=GORS-sparse(1:ng,1:ng,sum(A2GORS,1)+sum(D2GORS,1),ng,ng);   
   DORS1=DORS-sparse(1:nd,1:nd,sum(A2DORS,1)'+sum(D2BORS,2)+sum(D2CORS,2)+sum(D2GORS,2),nd,nd);
   
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
   
   AMORS=[AORS1,   A2CORS, A2GORS, A2DORS, A2BORS;
       A2CORS',  CORS1,  M2FR.C2G, D2CORS', M2FR.C2B;
       A2GORS', M2FR.C2G', GORS1,  D2GORS', M2FR.G2B;
       A2DORS', D2CORS, D2GORS,  DORS1,  D2BORS;
       A2BORS', M2FR.C2B',M2FR.G2B',D2BORS', zeros(nb,nb)];
   
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
   
    Aw = sparse(1:Nsum,1:Nsum,CMP.Cw,Nsum,Nsum)*AMW;
    Aoil = sparse(1:Nsum,1:Nsum,CMP.Co,Nsum,Nsum)*AMO;
    Aoilrs = sparse(1:Nsum,1:Nsum,CMP.Cor,Nsum,Nsum)*AMORS;
    Agas = sparse(1:Nsum,1:Nsum,CMP.Cg,Nsum,Nsum)*AMG;        
            
   [dItime,dItimeW,dItimeO,dItimeORs] = NR_TimeBO(GEOM.dV,Sw,So,Sw0,So0,CMP,dt);         
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
      
    AM = AM + Aw + Aoil + Aoilrs + Agas; 
    
    Qz1 = sparse(WELL.Won(Qf~=0,3),ones(1,size(WELL.Won(Qf~=0),1)),Qf(Qf~=0),nw,1);
    Qm1 = sparse(WELL.Won(Qf~=0,3),ones(1,size(WELL.Won(Qf~=0),1)),QQ.QQm(WELL.Won(Qf~=0),1),nw,1);
    Qc1 = sparse(WELL.WonC(Qf(WELL.WonC(:,3))~=0,3),ones(1,size(WELL.WonC(Qf(WELL.WonC(:,3))~=0,3),1)),QQ.QQc(WELL.WonC(Qf(WELL.WonC(:,3))~=0,1)),nw,1);
    Qg1 = sparse(WELL.WonG(Qf(WELL.WonG(:,3))~=0,3),ones(1,size(WELL.WonG(Qf(WELL.WonG(:,3))~=0,3),1)),QQ.QQg(WELL.WonG(Qf(WELL.WonG(:,3))~=0,1)),nw,1);
    Qd1 = sparse(WELL.WonD(Qf(WELL.WonD(:,3))~=0,3),ones(1,size(WELL.WonD(Qf(WELL.WonD(:,3))~=0,3),1)),QQ.QQd(WELL.WonD(Qf(WELL.WonD(:,3))~=0,1)),nw,1);
   
    Q1 = -Qz1+Qm1+Qc1+Qg1+Qd1;
    dPt = [BM',Q1(Qf~=0)']/[AM,W2M2;W2M1,W2M3];  %/////////////
    dPt = dPt';
      
    dPw(Qf~=0) = dPt(Nsum+1:end);
   
    [QQwoBND] = QIterGYBO(dPt,QQwoBND,WBND,BXYZ);
    
    [QQ.QQmwog] = QIterBO(QQ.QQmwog,W6,Wo,Wg,dPt(VEC.va),WELL.Won(:,1),dPw(WELL.Won(:,3)),na); 
    [QQ.QQcwog] = QIterBO(QQ.QQcwog,W6C,WoC,WgC,dPt(VEC.vc),WELL.WonC(:,1),dPw(WELL.WonC(:,3)),nc); 
    [QQ.QQgwog] = QIterBO(QQ.QQgwog,W6G,WoG,WgG,dPt(VEC.vg),WELL.WonG(:,1),dPw(WELL.WonG(:,3)),ng); 
    [QQ.QQdwog] = QIterBO(QQ.QQdwog,W6D,WoD,WgD,dPt(VEC.vd),WELL.WonD(:,1),dPw(WELL.WonD(:,3)),nd); 
   
 %  if isempty(RC.Cr)~=0 || isempty(RC.Gr)~=0
 %    flag_gim = 0;
 %  else
    Phi(1:Nsum,1) = Phi(1:Nsum,1) + dPt(1:Nsum); 
    Phi(1:Nsum,2) = Phi(1:Nsum,2) + dPt(1:Nsum); 
    Phi(1:Nsum,3) = Phi(1:Nsum,3) + dPt(1:Nsum);
    Pi(1:Nsum,1) = Pi(1:Nsum,1) + dPt(1:Nsum); 
    Pw = Pw + dPw;
 
    Fwater = AMW*Phi(1:Nsum,1) - dItimeW + [QQ.QQmwog(:,1);QQ.QQcwog(:,1);QQ.QQgwog(:,1);QQ.QQdwog(:,1);zeros(nb,1)] + [QQwoBND.Qmw;zeros(nc,1);zeros(ng,1);QQwoBND.Qdw;QQBND.Qb - WBND.b1gb.*dPt(VEC.vb)];
    FoilRs = AMORS*Phi(1:Nsum,2) - dItimeORs + CMP.Rs(:,2).*([QQ.QQmwog(:,2);QQ.QQcwog(:,2);QQ.QQgwog(:,2);QQ.QQdwog(:,2);zeros(nb,1)] + [QQwoBND.Qmo;zeros(nc,1);zeros(ng,1);QQwoBND.Qdo;QQBND.Qb - WBND.b1gb.*dPt(VEC.vb)]);
    Foil = AMO*Phi(1:Nsum,2) - dItimeO + [QQ.QQmwog(:,2);QQ.QQcwog(:,2);QQ.QQgwog(:,2);QQ.QQdwog(:,2);zeros(nb,1)] + [QQwoBND.Qmo;zeros(nc,1);zeros(ng,1);QQwoBND.Qdo;zeros(nb,1)];
     
    Sw(VSAT.vg,1) = Swl(VSAT.vg) +  (Fwater(VSAT.vg) - SGM.Cwp(VSAT.vg).*dPt(VSAT.vg))./SGM.Cwsw(VSAT.vg);
    So(VSAT.vg,1) = Sol(VSAT.vg) + (Foil(VSAT.vg) - SGM.Cop(VSAT.vg).*dPt(VSAT.vg))./SGM.Coso(VSAT.vg);
    Pb(VSAT.vg) = Pi(VSAT.vg);
    
    Sw(VSAT.vp) = Swl(VSAT.vp) + (Foil(VSAT.vp) - SGM.Cop(VSAT.vp).*dPt(VSAT.vp)).*SGM.Cgpb(VSAT.vp) - (FoilRs(VSAT.vp) - SGM.Cgp(VSAT.vp).*dPt(VSAT.vp)).*SGM.Copb(VSAT.vp);
    So(VSAT.vp) = 1 - Sw(VSAT.vp);
    Pb(VSAT.vp) = Pbl(VSAT.vp) + (FoilRs(VSAT.vp) - SGM.Cgp(VSAT.vp).*dPt(VSAT.vp)).*SGM.Cosw(VSAT.vp) - (Foil(VSAT.vp) - SGM.Cop(VSAT.vp).*dPt(VSAT.vp)).*SGM.Cgsw(VSAT.vp);
  
    vg1 = Pb > Pi + 1.e-10;
    vg = Pb(VSAT.vp) > Pi(VSAT.vp) + 1.e-10;     
    vp1 = So + Sw > 1 + 1.e-3;   
    vp = So(VSAT.vg) + Sw(VSAT.vg) > 1 + 1.e-3;  
    VSAT.vp(vg) = []; 
    VSAT.vg(vp) = []; 
    VSAT.vg = [VSAT.vg,VEC.v(vg1)]; VSAT.vp = [VSAT.vp,VEC.v(vp1)]; 
    
    Pb(VEC.v(vg1)) = Pi(VEC.v(vg1));
    So(VEC.v(vp1)) = 1 - Sw(VEC.v(vp1));
    
   %dBoP = - CMP.Bo(VEC.v(vg1),2)*PR.zc(2)./(1 + PR.zc(2)*(Pi(VEC.v(vg1)) - Pb(VEC.v(vg1))));
   % dBgP = - CMP.Bg(VEC.v(vg1),2)*PR.zc(3)./(1 + PR.zc(3)*(Pi(VEC.v(vg1)) - Pb(VEC.v(vg1))));
   % Bg = CMP.Bg(VEC.v(vg1),2) + dBgP.*dPt(VEC.v(vg1));
  %  Bo = CMP.Bo(VEC.v(vg1),2) + dBoP.*dPt(VEC.v(vg1));
  %  Rs = PR.rs.*Pi(VEC.v(vg1));
        
  %  So(VEC.v(vg1)) = (CMP.Rs(VEC.v(vg1),2)./CMP.Bo(VEC.v(vg1),2).*Sol(VEC.v(vg1)) + (Sw(VEC.v(vg1))- 1)./Bg)./(Rs./Bo - 1./Bg);
    %So(VEC.v(vg1)) = Sol(VEC.v(vg1)).*ds;
    
    So(VEC.v(vg1)) = Sol(VEC.v(vg1)) - 0.00001;
    Pb(VEC.v(vp1)) = Pbl(VEC.v(vp1)) - 0.00001;  
  
  %  QQ.QQm = CMP.Cw(VEC.va).*QQ.QQmwog(:,1)+(CMP.Co(VEC.va)+CMP.Rs(VEC.va,2).*CMP.Cor(VEC.va)).*QQ.QQmwog(:,2)+CMP.Cg(VEC.va).*QQ.QQmwog(:,3);
  %  QQ.QQc = CMP.Cw(VEC.vc).*QQ.QQcwog(:,1)+(CMP.Co(VEC.vc)+CMP.Rs(VEC.vc,2).*CMP.Cor(VEC.vc)).*QQ.QQcwog(:,2)+CMP.Cg(VEC.vc).*QQ.QQcwog(:,3);
  %  QQ.QQg = CMP.Cw(VEC.vg).*QQ.QQgwog(:,1)+(CMP.Co(VEC.vg)+CMP.Rs(VEC.vg,2).*CMP.Cor(VEC.vg)).*QQ.QQgwog(:,2)+CMP.Cg(VEC.vg).*QQ.QQgwog(:,3);
  % QQ.QQd = CMP.Cw(VEC.vd).*QQ.QQdwog(:,1)+(CMP.Co(VEC.vd)+CMP.Rs(VEC.vd,2).*CMP.Cor(VEC.vd)).*QQ.QQdwog(:,2)+CMP.Cg(VEC.vd).*QQ.QQdwog(:,3);
    
    QQ.QQm = QQ.QQm + sparse(WELL.Won(:,1),ones(1,nw),W1.*(dPw(WELL.Won(:,3))-dPt(VEC.va(WELL.Won(:,1)))),na,1);
    QQ.QQc = QQ.QQc + sparse(WELL.WonC(:,1),ones(1,size(W1C,1)),W1C.*(dPw(WELL.WonC(:,3))-dPt(VEC.vc(WELL.WonC(:,1)))),nc,1);
    QQ.QQg = QQ.QQg + sparse(WELL.WonG(:,1),ones(1,size(W1G,1)),W1G.*(dPw(WELL.WonG(:,3))-dPt(VEC.vg(WELL.WonG(:,1)))),ng,1);
    QQ.QQd = QQ.QQd + sparse(WELL.WonD(:,1),ones(1,size(W1D,1)),W1D.*(dPw(WELL.WonD(:,3))-dPt(VEC.vd(WELL.WonD(:,1)))),nd,1);
    
    Qm2 = QBild(QQ.QQm(WELL.Won(:,1)),QQ.QQmwog(WELL.Won(:,1),:),WELL.Uf(WELL.Won(:,3),ft),WELL.Won(:,1),dt,WELL.Won(:,3),nw,W1);
    Qc2 = QBild(QQ.QQc(WELL.WonC(:,1)),QQ.QQcwog(WELL.WonC(:,1),:),WELL.Uf(WELL.WonC(:,3),ft),WELL.WonC(:,1),dt,WELL.WonC(:,3),nw,W1C);
    Qg2 = QBild(QQ.QQg(WELL.WonG(:,1)),QQ.QQgwog(WELL.WonG(:,1),:),WELL.Uf(WELL.WonG(:,3),ft),WELL.WonG(:,1),dt,WELL.WonG(:,3),nw,W1G);       
    Qd2 = QBild(QQ.QQd(WELL.WonD(:,1)),QQ.QQgwog(WELL.WonD(:,1),:),WELL.Uf(WELL.WonD(:,3),ft),WELL.WonD(:,1),dt,WELL.WonD(:,3),nw,W1D);    
    flag_gim=sum(abs(dPt(1:na+nc+ng+nd)./Pi(1:na+nc+ng+nd))>=1e-6)~=0;    
  % end;
  end; 
  TL = sparse(1:na,1:na,CMP.Cw(VEC.va),na,na)*TW + sparse(1:na,1:na,CMP.Co(VEC.va),na,na)*TO + sparse(1:na,1:na,CMP.Cor(VEC.va),na,na)*TORS + sparse(1:na,1:na,CMP.Cg(VEC.va),na,na)*TGG;        
  DL = sparse(1:nd,1:nd,CMP.Cw(VEC.vd),nd,nd)*DW + sparse(1:nd,1:nd,CMP.Co(VEC.vd),nd,nd)*DO;% + sparse(1:nd,1:nd,CMP.Cor(VEC.vd),nd,nd)*DORS + sparse(1:nd,1:nd,CMP.Cg(VEC.vd),nd,nd)*DG;
     
 
   
   
