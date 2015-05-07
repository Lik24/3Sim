function [MPi,DPi,MPhi,DPhi,MSw,DSw,Cp,CMP,sQm2,sQd2,QQ]=Sat_fast2(Phi,Pi,Cp,Sw,Sw0,RC,KWOG,KWOG_GY,WELL,BFRACM,...
    nw,Qzm1,dt,Qf,GEOM,CMP,BB,TRM,M2FR,PR,fp,BXYZ,QQBND,QQoBND,CR_cr,Qm,Qd,QQ,ft,GYData)

na=RC.na;
nc=RC.nc;
ng=RC.ng;
nd=RC.nd;
nb=RC.nb;
vad=RC.ADr;
b1wb=sparse(nb,1);
va=1:na;  VEC.va = va;
          VEC.vc = zeros(0,1);
          VEC.vg = zeros(0,1);
v_d = na+nc+ng+1:na+nc+ng+nd;
vd=na+1:na+nd; VEC.vd = vd;
v_b = na+nc+ng+nd+1:na+nc+ng+nd+nb;
vb=na+nd+1:na+nd+nb; VEC.vb = vb;

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
dPw = zeros(nw,1);
dPt = zeros(na+nd,1);
flag_gim=1;
kj = 0;
while flag_gim==1 && kj<10
     kj=kj+1;
     [SGM,CMP]=SGimMD2([GEOM.dV(va);GEOM.dV(vd)],Sw,PR.zc,Pi,dPt,dt,CMP,[va,v_d]);
     [WBND,QQBND]=GY_bild2(GYData,Pi([va,vd],1),dPt,KWOG_GY,Cp([va,vd],1),RC,TRM.Txyz_GY_A,TRM.Txyz_GY_D,PR,CMP,VEC,QQBND);
     [W1,W6,Wo,W7]=Well_MKT2(dPt(va),dPw(WELL.Won(:,3)),WELL.Won,WELL.Uf(WELL.Won(:,3),ft+1),Cp(va,1),PR.aw,PR.as,PR,WELL.CpW(WELL.Won(:,3),ft+1),CMP,KWOG,va);% W1 - проводимость по всей жидкости, W6 - только дл€ воды, W7 - полимер 
     [W1D,W6D,WoD,W7D]=Well_MKT2(dPt(vd),dPw(WELL.WonD(:,3)),WELL.WonD,WELL.Uf(WELL.WonD(:,3),ft+1),Cp(vd,1),PR.tw,PR.ts,PR,WELL.CpW(WELL.WonD(:,3),ft+1),CMP,KWOG,v_d);

     b1wm=sparse(WELL.Won(:,1),ones(1,size(WELL.Won,1)),-W1.*dPw(WELL.Won(:,3)),na,1);
     b1wd=sparse(WELL.WonD(:,1),ones(1,size(WELL.WonD,1)),-W1D.*dPw(WELL.WonD(:,3)),nd,1);
     
     W2M=sparse(WELL.Won(:,3),WELL.Won(:,1),W1,nw,na);
     W2D=sparse(WELL.WonD(:,3),WELL.WonD(:,1),W1D,nw,nd);
     W2B=sparse(nw,nb);
     %Qf=Qzm1<0; ????????
     b1wm=b1wm.*(sum(W2M(Qf==0,:),1)~=0)';
     b1wd=b1wd.*(sum(W2D(Qf==0,:),1)~=0)';
     
     QQ.QQm = QQ.QQm + sparse(WELL.Won(:,1),ones(1,nw),W1.*(dPw(WELL.Won(:,3))-dPt(va(WELL.Won(:,1)))),na,1);
     QQ.QQd = QQ.QQd + sparse(WELL.WonD(:,1),ones(1,size(W1D,1)),W1D.*(dPw(WELL.WonD(:,3))-dPt(vd(WELL.WonD(:,1)))),nd,1);
    
     b1wm = b1wm - QQ.QQm;
     b1wd = b1wd - QQ.QQd; 
     
     WM1=[W2M,W2D,W2B];
     WM2=WM1';
     W3vec=sparse(WELL.Won(:,3),1,W1,nw,1)+sparse(WELL.WonD(:,3),1,W1D,nw,1);
     WM3=-sparse(1:nw,1:nw,W3vec,nw,nw);
     
     [TW,TO,TP]=Potok_MKT2(TRM.TTM,Phi(va,:),KWOG,Cp(va,1),PR,RC.Arc2,fp,PR.kms(1),GEOM.L,CMP);
     [DW,DO,DP]=Potok_Tube2(TRM.TD,Phi(vd,:),KWOG,Cp(vd,1),PR,fp,PR.kms(4),GEOM.Ld,RC.Dr2,RC.Dc2,nd,CMP,v_d);

     [A2DW,A2DO,A2DP]=Obmen_T2M2(M2FR.A2D,Phi(va,1),Phi(vd,1),GEOM.K(:,1),PR,Cp(va,1),Cp(vd,1),KWOG,CMP,va,v_d);
     [A2BW,A2BO,A2BP]=Obmen_T2M2(M2FR.A2B,Phi(va,1),Phi(vb,1),ones(na,1),PR,Cp(va,1),Cp(vb,1),KWOG,CMP,va,vb);
     [D2BW,D2BO,D2BP]=Obmen_T2M2(M2FR.D2B,Phi(vd,1),Phi(vb,1),GEOM.K(:,1),PR,Cp(vd,1),Cp(vb,1),KWOG,CMP,v_d,vb);
  
     A1 = -sparse(WELL.Won(:,1),WELL.Won(:,1),W1,na,na)-sparse(1:na,1:na,SGM.Clp(va)+sum(WBND.b1gm(:,1:2),2),na,na);
     D1 = -sparse(WELL.WonD(:,1),WELL.WonD(:,1),W1D,nd,nd)-sparse(1:nd,1:nd,SGM.Clp(vd)+sum(WBND.b1gd(:,1:2),2),nd,nd);                       %ћатрица коэф. дл€ двойной пор.
     B1 = -sparse(1:nb,1:nb,sum(A2BW,1)+sum(D2BW,1)+SGM.Clp(vb)'+WBND.b1gb',nb,nb);
        
   AW1=TW-sparse(1:na,1:na,sum(A2BW,2)+sum(A2DW,2),na,na);  
   DW1=DW-sparse(1:nd,1:nd,sum(A2DW,1)'+sum(D2BW,2),nd,nd);
      
   AO1=TO-sparse(1:na,1:na,sum(A2BO,2)+sum(A2DO,2),na,na);  
   DO1=DO-sparse(1:nd,1:nd,sum(A2DO,1)'+sum(D2BO,2),nd,nd);
    
   AMW=[AW1,   A2DW, A2BW;
        A2DW', DW1,  D2BW;
        A2BW', D2BW', B1];
       
   AMO=[AO1,   A2DO, A2BO;
        A2DO', DO1,  D2BO;
        A2BO', D2BO',zeros(nb,nb)];
    
     [dItime,dItimeW] = NR_TimeMD2([GEOM.dV(va);GEOM.dV(vd)],Sw,Sw0,CMP,dt,[va,v_d]);
     Awater = sparse(1:na+nd+nb,1:na+nd+nb,CMP.Cw([va,v_d]),na+nd+nb,na+nd+nb)*AMW;
     
     dPhi = Awater*Phi(:,1) + AMO*Phi(:,2);   
 
     bl = [b1wm;b1wd;b1wb] - [QQBND.Qm;QQBND.Qd;QQBND.Qb] + dItime - dPhi;
     
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
     bld=bl(vd)-Blcd1/dt-Blgd1/dt;
     Bl=[blm;bld;bl(vb)];
    
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
   
     Qm1 = sparse(WELL.Won(Qf~=0,3),ones(1,size(WELL.Won(Qf~=0),1)),QQ.QQm(Qf~=0),nw,1); 
     Qd1 = sparse(WELL.WonD(Qf(WELL.WonD(:,3))~=0,3),ones(1,size(WELL.WonD(Qf(WELL.WonD(:,3))~=0,3),1)),QQ.QQd(Qf(WELL.WonD(:,3))~=0),nw,1);
     
     AM = [A1,  zeros(na,nd), zeros(na,nb);
           zeros(nd,na), D1, zeros(nd,nb);
           zeros(nb,na), zeros(nb,nd), B1];
           
     AM = AM + Awater + AMO;
     
     dPt=[Bl',-Qzm1(Qf~=0)'+Qm1(Qf~=0)'+Qd1(Qf~=0)']/[AM,WM2;WM1,WM3];
     dPt = dPt';
     dPw(Qf~=0) = dPt(na+nd+nb+1:end);
     
   [QQoBND] = QIterGY2(dPt,QQoBND,WBND,BXYZ);
    
   [QQ.QQmwo]=QIter2(QQ.QQmwo,W6,Wo,dPt(va),WELL.Won(:,1),dPw(WELL.Won(:,3)),na);
   [QQ.QQdwo]=QIter2(QQ.QQdwo,W6D,WoD,dPt(vd),WELL.WonD(:,1),dPw(WELL.WonD(:,3)),nd); 
      
   Phi(:,1) = Phi(:,1) + dPt;   %//////////////
   Phi(:,2) = Phi(:,2) + dPt;   %//////////////
 
   Bw=[Bwc1;Bwcd1]/dt+[Bwg1;Bwgd1]/dt;
   Bo=[Boc1;Bocd1]/dt+[Bog1;Bogd1]/dt;
     
   Fwater = AMW*Phi(1:na+nd+nb,1) - dItimeW + [QQ.QQmwo(:,1);QQ.QQdwo(:,1)] + [QQoBND.Qmw;QQoBND.Qdw] + Bw;
   
   Sw = Sw + (Fwater - SGM.Cwp.*dPt )./SGM.Cwsw;
      
   Pi = Pi + dPt;
   flag_gim=sum(abs(dPt./Pi)>=1e-6)~=0;
end
     MSw=Sw(va); MPi=Pi(va); MPhi=Phi(va,:);
     DSw=Sw(vd); DPi=Pi(vd); DPhi=Phi(vd,:);
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
 %   Sw=Sw.*(Sw>=0).*(Sw<=1)+(Sw>1); 
 %   Cp=Cp.*(Cp>=0).*(Cp<=1)+(Cp>1);
    sQm2 = QBild(QQ.QQm(WELL.Won(:,3)),QQ.QQmwo(WELL.Won(:,3),:),WELL.Uf(WELL.Won(:,3),ft+1),WELL.Won(:,1),dt,WELL.Won(:,3),nw,W1);
    sQd2 = QBild(QQ.QQd(WELL.WonD(:,3)),QQ.QQdwo(WELL.WonD(:,3),:),WELL.Uf(WELL.WonD(:,3),ft+1),WELL.WonD(:,1),dt,WELL.WonD(:,3),nw,W1D);
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