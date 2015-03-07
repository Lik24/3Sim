function [Sw,So,Cp,sQm2,sQd2,QWOG]=Sat_fast_3(Phi,Pi,Cp,Sw,So,Sw0,So0,RC,KWOG,WELL,BFRACM,...
    nw,Qzm1,dt,Qf,GEOM,SGMMD,BB,TRM,M2FR,PR,fp,BOUNDFL,QBOUND,CR_cr,Qm,Qd,QWOG,ft)

na=RC.na;
nc=RC.nc;
ng=RC.ng;
nd=RC.nd;
nb=RC.nb;
vad=RC.ADr;
b1wb=sparse(nb,1);
dPw = zeros(size(Qzm1(:,1),1),1);

va=1:na;  VEC.va = va;
vc=na+1:na+nc; VEC.vc = zeros(0,1);
vg=na+nc+1:na+nc+ng; VEC.vg = zeros(0,1);
vd=na+nc+ng+1:na+nc+ng+nd; VEC.vd = vd;
vb=na+nc+ng+nd+1:na+nc+ng+nd+nb; VEC.vb = vb;

v1a=zeros(na,1);
v1a(RC.ACr)=1;
r1a=find(v1a==1); %€чейки матрицы с трещинами
   
v2a=zeros(na,1);
v2a(RC.AGr)=1;
r2a=find(v2a==1); %€чейки матрицы с горизонтальными трещинами

v1d=zeros(nd,1);
v1d(RC.DCr)=1;
r1d=find(v1d==1);
    
v2d=zeros(nd,1);
v2d(RC.DGr)=1;
r2d=find(v2d==1);
flag_gim=1;
kj = 0;
while flag_gim==1 && kj<10
     kj=kj+1;
     
     [W1,W6,Wo,Wg,W7]=Well_MKT(WELL.Won,WELL.Uf(WELL.WonM,ft+1),Sw(va,1),Cp(va,1),PR.aw,PR.as,PR,WELL.CpW(WELL.WonM,ft+1),SGMMD);% W1 - проводимость по всей жидкости, W6 - только дл€ воды, W7 - полимер 
     [W1D,W6D,WoD,WgD,W7D]=Well_MKT(WELL.WonD,WELL.Uf(WELL.WonD(:,3),ft+1),Sw(vd,1),Cp(vd,1),PR.tw,PR.ts,PR,WELL.CpW(WELL.WonD(:,3),ft+1),SGMMD);
  
     b1wm=sparse(WELL.Won(:,1),ones(1,size(WELL.Won,1)),-W1.*dPw(WELL.Won(:,3)),na,1);
     b1wd=sparse(WELL.WonD(:,1),ones(1,size(WELL.WonD,1)),-W1D.*dPw(WELL.WonD(:,3)),nd,1);
     
     W2M=sparse(WELL.Won(:,3),WELL.Won(:,1),W1,nw,na);
     W2D=sparse(WELL.WonD(:,3),WELL.WonD(:,1),W1D,nw,nd);
     W2B=sparse(nw,nb);
     %Qf=Qzm1<0; ????????
     b1wm=b1wm.*(sum(W2M(Qf==0,:),1)~=0)';
     b1wd=b1wd.*(sum(W2D(Qf==0,:),1)~=0)';
     
     b1wm = b1wm - sparse(WELL.Won(:,1),ones(1,size(WELL.Won,1)),QWOG.QQm,na,1);
     b1wd = b1wd - sparse(WELL.WonD(:,1),ones(1,size(WELL.WonD,1)),QWOG.QQd,nd,1); 
     
     WM1=[W2M,W2D,W2B];
     WM2=WM1';
     W3vec=sparse(WELL.Won(:,3),1,W1,nw,1)+sparse(WELL.WonD(:,3),1,W1D,nw,1);
     WM3=-sparse(1:nw,1:nw,W3vec,nw,nw);
     
     [TL,TW,TO,TGG,TWCgs,TOCgs,TP]=Potok_MKT(TRM.TTM,Phi(va,:),KWOG,Cp(va,1),PR,RC.Arc2,fp,PR.kms(1),GEOM.L,SGMMD);
     [DL,DW,DO,DG,DWCgs,DOCgs,DP]=Potok_Tube(TRM.TD,Phi(vd,:),Sw(vd,1),Cp(vd,1),PR,fp,PR.kms(4),GEOM.Ld,RC.Dr2,RC.Dc2,nd,SGMMD);
  
     [A2DL,A2DW,A2DO,A2DG,A2DWCgs,A2DOCgs,A2DP]=Obmen_T2M(M2FR.A2D,Pi(va,1),Pi(vd,1),Sw(va,1),Sw(vd,1),KWOG.K(:,1),PR,Cp(va,1),Cp(vd,1),SGMMD);
     [A2BL,A2BW,A2BO,A2BG,A2BWCgs,A2BOCgs,A2BP]=Obmen_T2M(M2FR.A2B,Pi(va,1),Pi(vb,1),Sw(va,1),Sw(vb,1),ones(na,1),PR,Cp(va,1),Cp(vb,1),SGMMD);
     [D2BL,D2BW,D2BO,D2BG,D2BWCgs,D2BOCgs,D2BP]=Obmen_T2M(M2FR.D2B,Pi(vd,1),Pi(vb,1),Sw(vd,1),Sw(vb,1),KWOG.K(:,1),PR,Cp(vd,1),Cp(vb,1),SGMMD);
  
     A1 = -sparse(WELL.Won(:,1),WELL.Won(:,1),W1,na,na)-sparse(1:na,1:na,SGMMD.Clp(va)+sum(BOUNDFL.b1gm(:,1:2),2),na,na);
     D1 = -sparse(WELL.WonD(:,1),WELL.WonD(:,1),W1D,nd,nd)-sparse(1:nd,1:nd,SGMMD.Clp(vd)+sum(BOUNDFL.b1gd(:,1:2),2),nd,nd);                       %ћатрица коэф. дл€ двойной пор.
     B1 = -sparse(1:nb,1:nb,sum(A2BL,1)+sum(D2BL,1)+SGMMD.Clp(vb)'+BOUNDFL.b1gb',nb,nb);
        
   AG1=TGG-sparse(1:na,1:na,sum(A2BG,2)+sum(A2DG,2),na,na); 
   DG1=DG-sparse(1:nd,1:nd,sum(A2DG,1)'+sum(D2BG,2),nd,nd);
        
   AW1Cgs=TWCgs-sparse(1:na,1:na,sum(A2BWCgs,2)+sum(A2DWCgs,2),na,na);  
   DW1Cgs=DWCgs-sparse(1:nd,1:nd,sum(A2DWCgs,1)'+sum(D2BWCgs,2),nd,nd);
      
   AO1Cgs=TOCgs-sparse(1:na,1:na,sum(A2BOCgs,2)+sum(A2DOCgs,2),na,na);  
   DO1Cgs=DOCgs-sparse(1:nd,1:nd,sum(A2DOCgs,1)'+sum(D2BOCgs,2),nd,nd);
   
   AW1=TW-sparse(1:na,1:na,sum(A2BW,2)+sum(A2DW,2),na,na);  
   DW1=DW-sparse(1:nd,1:nd,sum(A2DW,1)'+sum(D2BW,2),nd,nd);
      
   AO1=TO-sparse(1:na,1:na,sum(A2BO,2)+sum(A2DO,2),na,na);  
   DO1=DO-sparse(1:nd,1:nd,sum(A2DO,1)'+sum(D2BO,2),nd,nd);
    
   AMG=[AG1,   A2DG, A2BG;
        A2DG', DG1,  D2BG;
        A2BG', D2BG',zeros(nb,nb)];
        
   AMW=[AW1,   A2DW, A2BW;
        A2DW', DW1,  D2BW;
        A2BW', D2BW', B1];
       
   AMO=[AO1,   A2DO, A2BO;
        A2DO', DO1,  D2BO;
        A2BO', D2BO',zeros(nb,nb)];

   AMWCgs=[AW1Cgs,   A2DWCgs,  A2BWCgs;
           A2DWCgs', DW1Cgs,   D2BWCgs;
           A2BWCgs', D2BWCgs', zeros(nb,nb)];
       
   AMOCgs=[AO1Cgs,   A2DOCgs, A2BOCgs;
           A2DOCgs', DO1Cgs,  D2BOCgs;
           A2BOCgs', D2BOCgs', zeros(nb,nb)];

    [dItime,dItimeW,dItimeO] = NR_Time(GEOM.dV,Sw,So,Sw0,So0,SGMMD,va,zeros(0,1),zeros(0,1),vd,vb,dt);
            
     dPhi = AMWCgs*Phi(:,1) + AMOCgs*Phi(:,2) + AMG*Phi(:,3);   
 
     bl = [b1wm;b1wd;b1wb] - [QBOUND.Qm;QBOUND.Qd;QBOUND.Qb] + dItime - dPhi;
     
     Blc1=zeros(na,1);      Blc1(RC.ACr)=BFRACM.Blc;
     Bwc1=zeros(na,1);      Bwc1(RC.ACr)=BFRACM.Bwc;
     Boc1=zeros(na,1);      Boc1(RC.ACr)=BFRACM.Boc;
     
     Blcd1=zeros(nd,1);     Blcd1(RC.DCr)=BFRACM.Blcd;
     Bwcd1=zeros(nd,1);     Bwcd1(RC.DCr)=BFRACM.Bwcd;
     Bocd1=zeros(nd,1);     Bocd1(RC.DCr)=BFRACM.Bocd;

     Blg1=zeros(na,1);      Blg1(RC.AGr)=BFRACM.Blg;
     Bwg1=zeros(na,1);      Bwg1(RC.AGr)=BFRACM.Bwg;
     Bog1=zeros(na,1);      Bog1(RC.AGr)=BFRACM.Bog;
     
     Blgd1=zeros(nd,1);     Blgd1(RC.DGr)=BFRACM.Blgd;
     Bwgd1=zeros(nd,1);     Bwgd1(RC.DGr)=BFRACM.Bwgd;
     Bogd1=zeros(nd,1);     Bogd1(RC.DGr)=BFRACM.Bogd;
     
     blm=bl(va)-Blc1/dt-Blg1/dt;%sparse(r1a,ones(sum(v1a),1),Blc(),na,1)
     bld=bl(vd-nc-ng)-Blcd1/dt-Blgd1/dt;
     Bl=[blm;bld;bl(vb-na-nd)];
    
     WM1=WM1(Qf~=0,:);
     WM2=WM2(:,Qf~=0);
     WM3=WM3(Qf~=0,Qf~=0);
     
     Qzm=zeros(size(Qzm1));
     Qzd=zeros(size(Qzm1));
      
     Qzm1(CR_cr(1,1).won(:,3))=0;
     Qzm1(CR_cr(1,2).won(:,3))=0;
     
     Qzm(CR_cr(1,1).won(:,3))=(Qm(CR_cr(1,1).won(:,3),1)+Qm(CR_cr(1,1).won(:,3),2))/dt;% - -;%?????????
     Qzd(CR_cr(1,2).won(:,3))=(Qd(CR_cr(1,2).won(:,3),1)+Qd(CR_cr(1,2).won(:,3),2))/dt;
     Qzm1=Qzm1+Qzm+Qzd;
     
     Qm1 = sparse(WELL.WonM(Qf~=0),ones(1,size(WELL.WonM(Qf~=0),1)),QWOG.QQm(Qf~=0),nw,1);
     Qd1 = sparse(WELL.WonD(Qf(WELL.WonD(:,3))~=0,3),ones(1,size(WELL.WonD(Qf(WELL.WonD(:,3))~=0,3),1)),QWOG.QQd(Qf(WELL.WonD(:,3))~=0),nw,1);
   
     AM = [A1,  zeros(na,nd), zeros(na,nb);
           zeros(nd,na), D1, zeros(nd,nb);
           zeros(nb,na), zeros(nb,nd), B1];
           
     AM = AM + AMWCgs + AMOCgs + AMG;
     
     dPt=[Bl',-Qzm1(Qf~=0)'+Qm1(Qf~=0)'+Qd1(Qf~=0)']/[AM,WM2;WM1,WM3];
     dPt = dPt';
     dPw(Qf~=0) = dPt(na+nd+nb+1:end);
   
   [QWOG.QQm,QQmwo]=QIter(QWOG.QQm,W6,Wo,Wg,dPt(va),WELL.Won(:,1),dPw(WELL.Won(:,3)),SGMMD);
   [QWOG.QQd,QQdwo]=QIter(QWOG.QQd,W6D,WoD,WgD,dPt(vd),WELL.WonD(:,1),dPw(WELL.WonD(:,3)),SGMMD); 
   
   QWOG.QQmwo(:,1) = QWOG.QQmwo(:,1) + sparse(WELL.Won(:,1),ones(1,size(WELL.Won,1)),QQmwo(:,1),na,1);
   QWOG.QQdwo(:,1) = QWOG.QQdwo(:,1) + sparse(WELL.WonD(:,1),ones(1,size(WELL.WonD,1)),QQdwo(:,1),nd,1);
   
   QWOG.QQmwo(:,2) = QWOG.QQmwo(:,2) + sparse(WELL.Won(:,1),ones(1,size(WELL.Won,1)),QQmwo(:,2),na,1);
   QWOG.QQdwo(:,2) = QWOG.QQdwo(:,2) + sparse(WELL.WonD(:,1),ones(1,size(WELL.WonD,1)),QQdwo(:,2),nd,1);
   
   Phi(1:na+nd,1) = Phi(1:na+nd,1) + dPt(1:na+nd);   %//////////////
   Phi(1:na+nd,2) = Phi(1:na+nd,2) + dPt(1:na+nd);   %//////////////
  %Phi(:,3) = Phi(:,3) + dPt(1:na+nc+ng+nd);
   
   Bw=[Bwc1;Bwcd1]/dt+[Bwg1;Bwgd1]/dt;
   Bo=[Boc1;Bocd1]/dt+[Bog1;Bogd1]/dt;
     
   Fwater = AMW*Phi(1:na+nd,1) - dItimeW + [QWOG.QQmwo(:,1);QWOG.QQdwo(:,1)] + Bw;
   Foil = AMO*Phi(1:na+nd,2) - dItimeO + [QWOG.QQmwo(:,2);QWOG.QQdwo(:,2)] + Bo;
   Sw([va,vd]) = Sw([va,vd]) + (Fwater - SGMMD.Cwp.*dPt(1:na+nd))./SGMMD.Cws;
   So([va,vd]) = 1 - Sw([va,vd]);
   
   [SGMMD]=SGim(GEOM.dV,Sw,So,SGMMD.Mp(:,2),PR.zc,SGMMD.Rs(:,2),PR.rs,SGMMD.Bwog,dPt,Pi,1,va,zeros(0,1),zeros(0,1),vd,vb,dt);
  
   Pi(va)=Pi(va) + dPt(va);
   Pi(vd)=Pi(vd-na-nc-ng) + dPt(vd-na-nc-ng);
   Pi(vb)=Pi(vb-na-nc-ng-nd) + dPt(vb-na-nc-ng-nd);
   flag_gim=sum(abs(dPt([va,vd])./Pi([va,vd]))>=1e-6)~=0;
end
     
   %  b_A2B=Soed2B(A2BW,A2BP,Pi,na,nb,va,vb);     % —в€зь пор с граничной областью
   %  b_D2B=Soed2B(D2BW,D2BP,Pi,nd,nb,vd,vb);     % —в€зь трещин с граничной областью

    % bwm=sparse(Won,ones(1,size(Won,1)),-W6.*(Pi(Won)-PwNl(WonM)),na,1)...
     %    -b1gm(:,3).*(Pi(va)-GYData.GY_Pxy)-b1gm(:,4).*(Pi(va)-GYData.GY_Pz)-b_A2B(:,1);
     
    % bwd=sparse(WoD(:,1),ones(1,size(WoD,1)),-W6D.*(Pi(WoD(:,1))-PwNl(WoD(:,3))),nd,1)...
   %      -b1gd(:,3).*(Pi(vd)-GYData.GY_Pxy(vad))-b1gd(:,4).*(Pi(vd)-GYData.GY_Pz(vad))-b_D2B(:,1);
    
%     bpm=sparse(Won,ones(1,size(Won,1)),-W7.*(Pi(Won)-PwNl(WonM)),na,1)...
%         -b1gm(:,5).*(Pi(va)-GYData.GY_Pxy)-b1gm(:,6).*(Pi(va)-GYData.GY_Pz)-b_A2B(:,2);
     
%     bpd=sparse(WoD(:,1),ones(1,size(WoD,1)),-W7D.*(Pi(WoD(:,1))-PwNl(WoD(:,3))),nd,1)...
%         -b1gd(:,5).*(Pi(vd)-GYData.GY_Pxy(vad))-b1gd(:,6).*(Pi(vd)-GYData.GY_Pz(vad))-b_D2B(:,2);
     
    % bw=[bwm;bwd];
%     bp=[bpm;bpd];
     %blm=bl(va)-Blc1/dt-sparse(r2a,ones(sum(v2a),1),Blg,na,1)/dt;%sparse(r1a,ones(sum(v1a),1),Blc(),na,1)
     %bld=bl(vd-nc-ng)-Blcd1/dt-sparse(r2d,ones(sum(v2d),1),Blgd,nd,1)/dt;
    
     
     %Bw=bw+sparse([r1a;r1d+na],ones(sum(v1a)+sum(v1d),1),[Bwc;Bwcd],na+nd,1)/dt+sparse([r2a;r2d+na],ones(sum(v2a)+sum(v2d),1),[Bwg;Bwgd],na+nd,1)/dt...
     %    -Cwp([va,vd]).*(Pi([va,vd])-Pi0([va,vd]))+Grw([va,vd]);
 %    Bp=bp+sparse([r1a;r1d+na],ones(sum(v1a)+sum(v1d),1),[Bpc;Bpcd],na+nd,1)/dt+sparse([r2a;r2d+na],ones(sum(v2a)+sum(v2d),1),[Bpg;Bpgd],na+nd,1)/dt...
 %        -Cp([va,vd]).*Cwp([va,vd]).*(Pi([va,vd])-Pi0([va,vd]))+Cp([va,vd]).*Grw([va,vd]);
     %tmp=sum(Bw);    

     %AM2=[TW-sparse(1:na,1:na,sum(A2DW,2)+Cwp(va)),A2DW;
     %     A2DW',DW-sparse(1:nd,1:nd,sum(A2DW,1)+Cwp(vd)')];
 %    AM3=[TP,A2DP;A2DP',DP];
     
     %Sw_old=Sw([va,vd]);
     %Sw([va,vd])=Sw([va,vd])+dt*(AM2*Pi([va,vd])+Bw)./Cws([va,vd]);

 %    Cp([va,vd])=Sw0.*Cp([va,vd])+dt*(AM3*Pi([va,vd])+Bp)./Cws([va,vd]);
     
 %    vad=[va,vd];
 %    v0=Sw(vad)==0;
 %    Cp(vad(v0==0))=Cp(vad(v0==0))./Sw(vad(v0==0));

    % Sw=Sw.*(Sw>=PR.aw(4)).*(Sw<=1-PR.aw(5))+(1-PR.aw(5)).*(Sw>1-PR.aw(5))+PR.aw(4).*(Sw<PR.aw(4));
    Sw=Sw.*(Sw>=0).*(Sw<=1)+(Sw>1); 
 %   Cp=Cp.*(Cp>=0).*(Cp<=1)+(Cp>1);
    
  sQm2 = QBild(QWOG.QQm(WELL.Won(:,3)),QWOG.QQmwo(WELL.Won(:,1),:),WELL.Uf(WELL.Won(:,3),ft+1),WELL.Won(:,1),dt,WELL.Won(:,3),nw,W1);
  sQd2 = QBild(QWOG.QQd(WELL.WonD(:,3)),QWOG.QQdwo(WELL.WonD(:,1),:),WELL.Uf(WELL.WonD(:,3),ft+1),WELL.WonD(:,1),dt,WELL.WonD(:,3),nw,W1D);
end

function b_A2B=Soed2B(A2BW,A2BP,Pi,na,nb,va,vb)
[r1,c1]=find(A2BW);
dP=(Pi(va(r1))-Pi(vb(c1)));
A2B=A2BW(r1+(c1-1)*na).*dP;
A2B1=sparse(r1,c1,A2B,na,nb);
b_A2B(:,1)=sum(A2B1,2);

A2B=A2BP(r1+(c1-1)*na).*dP;
A2B1=sparse(r1,c1,A2B,na,nb);
b_A2B(:,2)=sum(A2B1,2);
end