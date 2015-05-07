function [BFRACM,Pi,Phi,SW,SW0,Pw,CMP,i,QQ]=fun2(RC,Pi,Phi,SW,Cp,PR,TRM,M2FR,WELL,Pw,dt,CR_rc,Qz,Qf,ndt,GEOM,Lc,CMP,ft,QQ,KWOG)
                                                                              
as=PR.as;
aw=PR.aw;
ts=PR.ts;
tw=PR.tw;
mu=PR.mu;
kms=PR.kms;
Ro=PR.Ro;
zc=PR.zc;
rs=PR.rs;
Na=RC.na;
nc=RC.nc;
ng=RC.ng;
Nd=RC.nd;
Uf = WELL.Uf(:,ft+1);
CpW = WELL.CpW(:,ft+1);
%C2G=sparse(nc,ng);     C2GL=C2G;  C2GW=C2G;

[v1,vc1,wom,r1,c1,r2,c2,rc_gy,rc_in_h,T_gy,T_in]=ext_cr(CR_rc(1,1));
[v2,vc2,wod,r1d,c1d,r2d,c2d,rc_gy_d,rc_in_hd,T_gy_d,TD_in]=ext_cr(CR_rc(1,2));
A2D=CR_rc(1,3).a2d;
r3=CR_rc(1,3).r;    c3=CR_rc(1,3).c;
dZ=CR_rc(1,4).dZ;

v_a=1:Na;
v_c=Na+1:Na+nc;
v_g=Na+nc+1:Na+nc+ng;
v_d=Na+nc+ng+1:Na+nc+ng+Nd;

na=sum(v1);
nd=sum(v2);

va=1:na;
vc=na+1:na+nc;
vg=na+nc+1:na+nc+ng;
vd=na+nc+ng+1:na+nc+ng+nd;
vb = zeros(1,0);
Pa=Pi(v_a);    Pc=Pi(v_c);    Pg=Pi(v_g);      Pd=Pi(v_d);
Phia=Phi(v_a,:); Phic=Phi(v_c,:); Phig=Phi(v_g,:);   Phid=Phi(v_d,:);

dVa=GEOM.dV(v_a);   dVc=GEOM.dV(v_c);   dVg=GEOM.dV(v_g);     dVd=GEOM.dV(v_d); Ke=GEOM.Ke;
MSw=SW(v_a);                               DSw=SW(v_d); 
MMp=CMP.Mp(v_a,:);                         DMp=CMP.Mp(v_d,:); 
MBw=CMP.Bw(v_a,:);                         DBw=CMP.Bw(v_d,:);
MBo=CMP.Bo(v_a,:);                         DBo=CMP.Bo(v_d,:);

dVa(v1~=1)=[];   dVd(v2~=1)=[];

Qz=Qz/dt;
nw=size(Qz,1);
Qf1 = Qf;
Qf=Qf(wom(:,3));

%%
vP1=Pa(RC.ACr)>=Pc(RC.ACc);
vP2=Pa(RC.AGr)>=Pg(RC.AGc);
    fl1=sum([SW(v_c);SW(v_g)])/size([SW(v_c);SW(v_g)],1)==1;
    fl2=sum([vP1;vP2])==0;
    Fl=fl1*fl2;
    
%%

Bwc=zeros(size(RC.ACr,1),1);   Boc=Bwc;   Blc=Bwc;
Bwg=zeros(size(RC.AGr,1),1);   Bog=Bwg;   Blg=Bwg;

Bwcd=zeros(size(RC.DCr,1),1);  Bocd=Bwcd; Blcd=Bwcd;
Bwgd=zeros(size(RC.DGr,1),1);  Bogd=Bwgd; Blgd=Bwgd;

sQ1=zeros(nw,5);
sQ2=zeros(nw,5);
sQm=zeros(nw,5);
sQd=zeros(nw,5);

Pj(:,1)=[Pa(v1==1);Pc;Pg;Pd(v2==1)];
Phj=[Phia(v1==1,:);Phic;Phig;Phid(v2==1,:)];

Sw=[MSw(v1==1);SW(v_c);SW(v_g);DSw(v2==1)];

CMPF.Bw=[MBw(v1==1,:);CMP.Bw(v_c,:);CMP.Bw(v_g,:);DBw(v2==1,:)];
CMPF.Bo=[MBo(v1==1,:);CMP.Bo(v_c,:);CMP.Bo(v_g,:);DBo(v2==1,:)];
CMPF.Mp=[MMp(v1==1,:);CMP.Mp(v_c,:);CMP.Mp(v_g,:);DMp(v2==1,:)];
CMPF.Mp0 = [CMP.Mp0(v1==1,:);CMP.Mp0(v_c,:);CMP.Mp0(v_g,:);CMP.Mp0(v2==1,:)];
CMPF.P0 = [CMP.P0(v1==1,:);CMP.P0(v_c,:);CMP.P0(v_g,:);CMP.P0(v2==1,:)];

Pgy=Pa(rc_gy(:,1));
Pgy2=Pd(rc_gy_d(:,1));
 
SW0=SW;

dPA = zeros(size(Pa));
dPD = zeros(size(Pd));

kfw=zeros(size(Sw));
kfo=zeros(size(Sw));

KfwM=KWOG.w(v_a); %water
KfoM=KWOG.o(v_a); %oil

KfwD=KWOG.w(v_d); %water
KfoD=KWOG.o(v_d); %oil

[bAl,ql,qo,qw]=Potok_GY20(T_gy,Pgy,Pa,rc_gy,KfwM,KfoM,v1,mu,Na,vc1,CMP);
[bDl,qld,qod,qwd]=Potok_GY20(T_gy_d,Pgy2,Pd,rc_gy_d,KfwD,KfoD,v2,mu,Nd,vc2,CMP);

QQm = QQ.QQm(v1==1); QQmwo = QQ.QQmwo(v1==1,:);
QQd = QQ.QQd(v2==1); QQdwo = QQ.QQdwo(v2==1,:);

wna=unique(wom(:,3));
wnc=unique(WELL.WonC(:,3));
wng=unique(WELL.WonG(:,3));
wnd=unique(wod(:,3));

fl=1;
i=0;
j_ndt=0;
fl2=0;
PjOld = Pj;
while i<ndt
    i = i + 1;
    dPj = zeros(size(Sw));
    dPw = zeros(size(Qz(:,1),1),1);
    kj = 0;
    flag_gim = 1;
    kfw(va)=Sat_cal(Sw(va),1-Sw(va),1,1,as,aw); %water
    kfo(va)=Sat_cal(Sw(va),1-Sw(va),2,1,as,aw); %oil
   
    kfw([vc,vg,vd])=Sat_cal(Sw([vc,vg,vd]),1-Sw([vc,vg,vd]),1,1,ts,tw); %water
    kfo([vc,vg,vd])=Sat_cal(Sw([vc,vg,vd]),1-Sw([vc,vg,vd]),2,1,ts,tw); %oil
    
    KfwM(v1==1)=kfw(va);      KfoM(v1==1)=kfo(va);
    KfwD(v2==1)=kfw(vd);      KfoD(v2==1)=kfo(vd);
    
    KFA=[kfw(va),kfo(va)];
    KFC=[kfw(vc),kfo(vc)];
    KFG=[kfw(vg),kfo(vg)];
    KFD=[kfw(vd),kfo(vd)];
    Sw0=Sw;

    CMPF.Bw(:,1) = CMPF.Bw(:,2);
    CMPF.Bo(:,1) = CMPF.Bo(:,2);
    CMPF.Mp(:,1) = CMPF.Mp(:,2);
     
    while flag_gim==1 && kj<10
        kj=kj+1;
        [SGMF,CMPF]=SGim2([dVa;dVc;dVg;dVd],Sw,zc,Pj,dPj,dt/ndt,CMPF);
        CMP.Bo([v_a(v1==1),v_c,v_g,v_d(v2==1)],2) = CMPF.Bo([va,vc,vg,vd],2);
        CMP.Bw([v_a(v1==1),v_c,v_g,v_d(v2==1)],2) = CMPF.Bw([va,vc,vg,vd],2);
        CMP.Mp([v_a(v1==1),v_c,v_g,v_d(v2==1)],2) = CMPF.Mp([va,vc,vg,vd],2);
        CMP.Cw([v_a(v1==1),v_c,v_g,v_d(v2==1)]) = CMPF.Cw([va,vc,vg,vd]);
        
        PjOld = Pj;
        [vPa1,vPc1,vPg1,vPd1,dPa,dPc,dPg,dPd]=pre_potok_2(Phi,Phj,RC,rc_in_h,rc_in_hd,Na,Nd,na,nd);

        [TW,TO,TP]=Potok_MKTF2(T_in,vPa1,Cp(v_a),mu,rc_in_h,Na,KfwM,KfoM,kms(1),dPa,GEOM.L,Ro,Ke,CMP,v1,v_a);
        [CW,CO]=Potok_TubeF2(TRM.TC,Pc,vPc1,kfw(vc),kfo(vc),Cp(v_c),PR,RC.Cr2,RC.Cc2,kms(2),dPc,Lc,nc,CMP,vc);
        [GW,GO]=Potok_TubeF2(TRM.TG,Pg,vPg1,kfw(vg),kfo(vg),Cp(v_g),PR,RC.Gr2,RC.Gc2,kms(3),dPg,GEOM.Lg,ng,CMP,vg);
        [DW,DO,DP]=Potok_MKTF2(TD_in,vPd1,Cp(v_d),mu,rc_in_hd,Nd,KfwD,KfoD,kms(4),dPd,GEOM.L,Ro,Ke,CMP,v2,v_d);
        
        Pa1=Pi(v_a(v1==1));
        Cp1=Cp(v_a(v1==1));
        
        Pd1=Pi(v_d(v2==1));
        Cp1d=Cp(v_d(v2==1));
       
        [A2CL,A2CW,A2CO]=Obmen_T2MF2(TRM.TA2C,Pa1,Pc,mu,Cp1,Cp(v_c),r1,RC.ACc,KFA,KFC,CMPF,va,vc);
        [A2GL,A2GW,A2GO]=Obmen_T2MF2(TRM.TA2G,Pa1,Pg,mu,Cp1,Cp(v_g),r2,RC.AGc,KFA,KFG,CMPF,va,vg);
        [~,A2DW,A2DO]=Obmen_T2MF2(A2D,Pa1,Pd1,mu,Cp1,Cp1d,r3,c3,KFA,KFD,CMPF,va,vd);
        
        [D2CL,D2CW,D2CO]=Obmen_T2MF2(TRM.TD2C,Pd1,Pc,mu,Cp1d,Cp(v_c),r1d,RC.DCc,KFD,KFC,CMPF,vd,vc);
        [D2GL,D2GW,D2GO]=Obmen_T2MF2(TRM.TD2G,Pd1,Pg,mu,Cp1d,Cp(v_g),r2d,RC.DGc,KFD,KFG,CMPF,vd,vg);
        
        [W1,W6,Wo,W7]=Well_MKTF2(dPj(va),dPw(wom(:,3)),wom(:,2),wom(:,1),Uf(wom(:,3)),Cp(v_a(v1==1)),mu,CpW(wom(:,3)),kfw(va),kfo(va),CMP,CMPF,QQ.QQm(wom(:,3)));
        [W1C,W6C,WoC,W7C]=Well_MKTF2(dPj(vc),dPw(WELL.WonC(:,3)),WELL.WonC(:,2),WELL.WonC(:,1),Uf(WELL.WonC(:,3)),Cp(v_c),mu,CpW(WELL.WonC(:,3)),kfw(vc),kfo(vc),CMP,CMPF,QQ.QQc);
        [W1G,W6G,WoG,W7G]=Well_MKTF2(dPj(vg),dPw(WELL.WonG(:,3)),WELL.WonG(:,2),WELL.WonG(:,1),Uf(WELL.WonG(:,3)),Cp(v_g),mu,CpW(WELL.WonG(:,3)),kfw(vg),kfo(vg),CMP,CMPF,QQ.QQg);
        [W1D,W6D,WoD,W7D]=Well_MKTF2(dPj(vd),dPw(wod(:,3)),wod(:,2),wod(:,1),Uf(wod(:,3)),Cp(v_d(v2==1)),mu,CpW(wod(:,3)),kfw(vd),kfo(vd),CMP,CMPF,QQ.QQd(wod(:,3)));
                
        AW1=TW-sparse(1:na,1:na,sum(A2CW,2)+sum(A2GW,2)+sum(A2DW,2),na,na);
        CW1=CW-sparse(1:nc,1:nc,sum(A2CW,1)+sum(D2CW,1),nc,nc);
        GW1=GW-sparse(1:ng,1:ng,sum(A2GW,1)+sum(D2GW,1),ng,ng);
        DW1=DW-sparse(1:nd,1:nd,sum(A2DW,1)'+sum(D2CW,2)+sum(D2GW,2),nd,nd);
        
        AO1=TO-sparse(1:na,1:na,sum(A2CO,2)+sum(A2GO,2)+sum(A2DO,2),na,na);
        CO1=CO-sparse(1:nc,1:nc,sum(A2CO,1)+sum(D2CO,1),nc,nc);
        GO1=GO-sparse(1:ng,1:ng,sum(A2GO,1)+sum(D2GO,1),ng,ng);
        DO1=DO-sparse(1:nd,1:nd,sum(A2DO,1)'+sum(D2CO,2)+sum(D2GO,2),nd,nd);
              
        AMW=[AW1,   A2CW, A2GW, A2DW;
            A2CW',  CW1,  M2FR.C2G, D2CW';
            A2GW', M2FR.C2G', GW1,  D2GW';
            A2DW', D2CW, D2GW,  DW1];
        
        AMO=[AO1,   A2CO, A2GO, A2DO;
            A2CO',  CO1,  M2FR.C2G, D2CO';
            A2GO', M2FR.C2G', GO1,  D2GO';
            A2DO', D2CO, D2GO,  DO1];
        
        A1 = - sparse(1:na,1:na,SGMF.Clp(va) + bAl,na,na)-sparse(wom(:,1),wom(:,1),W1,na,na);
        C1 = - sparse(1:nc,1:nc,SGMF.Clp(vc)',nc,nc)-sparse(WELL.WonC(:,1),WELL.WonC(:,1),W1C,nc,nc);
        G1 = - sparse(1:ng,1:ng,SGMF.Clp(vg)',ng,ng)-sparse(WELL.WonG(:,1),WELL.WonG(:,1),W1G,ng,ng);
        D1 = - sparse(1:nd,1:nd,SGMF.Clp(vd) + bDl,nd,nd)-sparse(wod(:,1),wod(:,1),W1D,nd,nd);
        
        AM = [A1, zeros(na,nc), zeros(na,ng), zeros(na,nd);
            zeros(nc,na), C1, zeros(nc,ng), zeros(nc,nd);
            zeros(ng,na), zeros(ng,nc), G1, zeros(ng,nd);
            zeros(nd,na), zeros(nd,nc), zeros(nd,ng), D1];
        
        Awater = sparse(1:na+nc+ng+nd,1:na+nc+ng+nd,CMPF.Cw,na+nc+ng+nd,na+nc+ng+nd)*AMW;
        AM = AM + Awater + AMO;
        
        W2M=sparse(wom(:,3),wom(:,1),W1,nw,na);
        W2C=sparse(WELL.WonC(:,3),WELL.WonC(:,1),W1C,nw,nc);
        W2G=sparse(WELL.WonG(:,3),WELL.WonG(:,1),W1G,nw,ng);
        W2D=sparse(wod(:,3),wod(:,1),W1D,nw,nd);
        
        WM1=[W2M,W2C,W2G,W2D];
        WM2=WM1';
        W3vec=sparse(wom(:,3),1,W1,nw,1)+sparse(WELL.WonC(:,3),1,W1C,nw,1)+sparse(WELL.WonG(:,3),1,W1G,nw,1)+sparse(wod(:,3),1,W1D,nw,1);
        WM3=-sparse(1:nw,1:nw,W3vec,nw,nw);
        %WM1(Qz~=0,:)
        % full([Qz(wn(Qf~=0)),W3vec(wn(Qf~=0)),Pw(wn(Qf~=0))])
        WM1=WM1(wom(Qf~=0,3),:);
        WM2=WM2(:,wom(Qf~=0,3));
        WM3=WM3(wom(Qf~=0,3),wom(Qf~=0,3));
        
        ba1=sparse(wom(:,1),ones(1,size(wom,1)),-W1.*dPw(wom(:,3)),na,1);
        bc1=sparse(WELL.WonC(:,1),ones(1,size(WELL.WonC,1)),-W1C.*dPw(WELL.WonC(:,3)),nc,1);
        bg1=sparse(WELL.WonG(:,1),ones(1,size(WELL.WonG,1)),-W1G.*dPw(WELL.WonG(:,3)),ng,1);
        bd1=sparse(wod(:,1),ones(1,size(wod,1)),-W1D.*dPw(wod(:,3)),nd,1);
        
        ba1=ba1.*(sum(W2M(Qz==0,:),1)~=0)';
        bc1=bc1.*(sum(W2C(Qz==0,:),1)~=0)';
        bg1=bg1.*(sum(W2G(Qz==0,:),1)~=0)';
        bd1=bd1.*(sum(W2D(Qz==0,:),1)~=0)';
  
        QQm = QQm + sparse(wom(:,1),ones(1,size(W1,1)),W1.*(dPw(wom(:,3))-dPj(va(wom(:,1)))),na,1);
        QQ.QQc = QQ.QQc + sparse(WELL.WonC(:,1),ones(1,size(W1C,1)),W1C.*(dPw(WELL.WonC(:,3))-dPj(vc(WELL.WonC(:,1)))),nc,1);
        QQ.QQg = QQ.QQg + sparse(WELL.WonG(:,1),ones(1,size(W1G,1)),W1G.*(dPw(WELL.WonG(:,3))-dPj(vg(WELL.WonG(:,1)))),ng,1);
        QQd = QQd + sparse(wod(:,1),ones(1,size(W1D,1)),W1D.*(dPw(wod(:,3))-dPj(vd(wod(:,1)))),nd,1);
        
        ba1 = ba1 - QQm;
        bc1 = bc1 - QQ.QQc;
        bg1 = bg1 - QQ.QQg;
        bd1 = bd1 - QQd;
        
        [dItime,dItimeW] = NR_Time2([dVa;dVc;dVg;dVd],Sw,Sw0,CMPF,dt/ndt);
        
        dPhj = Awater*Phj(:,1) + AMO*Phj(:,2);
        
        BC = [ba1 - ql;bc1;bg1;bd1 - qld] + dItime - dPhj;
       
        Qz1 = sparse(wom(Qf~=0,3),ones(1,size(wom(Qf~=0,3),1)),Qz(wom(Qf~=0,3)),nw,1);
        Qm1 = sparse(wom(Qf1(wom(:,3))~=0,3),ones(1,size(wom(Qf1(wom(:,3))~=0,3),1)),QQm(wom(Qf1(wom(:,3))~=0,1)),nw,1);
        Qc1 = sparse(WELL.WonC(Qf1(WELL.WonC(:,3))~=0,3),ones(1,size(WELL.WonC(Qf1(WELL.WonC(:,3))~=0,3),1)),QQ.QQc(WELL.WonC(Qf1(WELL.WonC(:,3))~=0,1)),nw,1);
        Qg1 = sparse(WELL.WonG(Qf1(WELL.WonG(:,3))~=0),ones(1,size(WELL.WonG(Qf1(WELL.WonG(:,3))~=0,3),1)),QQ.QQg(WELL.WonC(Qf1(WELL.WonG(:,3))~=0,1)),nw,1);
        Qd1 = sparse(wod(Qf1(wod(:,3))~=0,3),ones(1,size(wod(Qf1(wod(:,3))~=0,3),1)),QQd(wod(Qf1(wod(:,3))~=0,1)),nw,1);
        
        Q11 = -Qz1+Qm1+Qc1+Qg1+Qd1;

        dPt=[BC',Q11(wna(Qf~=0))']/[AM,WM2;WM1,WM3];
        dPt = dPt';
        dPj(:,1)=dPt(1:na+nc+ng+nd);
        dPA(v1==1)=dPj(va);
        dPD(v2==1)=dPj(vd);
        dPw(wom(Qf~=0,3))=dPt(na+nc+ng+nd+1:end);
        
        Phj(1:na+nc+ng+nd,1) = Phj(1:na+nc+ng+nd,1) + dPj;   %//////////////
        Phj(1:na+nc+ng+nd,2) = Phj(1:na+nc+ng+nd,2) + dPj;
              
        [bAl,ql,qw]=Potok_GY2(T_gy,Pgy,Pa,dPA,rc_gy,KfwM,KfoM,v1,mu,Na,vc1,CMP,qw,ql,va);
        [bDl,qld,qwd]=Potok_GY2(T_gy_d,Pgy2,Pd,dPD,rc_gy_d,KfwD,KfoD,v2,mu,Nd,vc2,CMP,qwd,qld,vd);
        
        [QQmwo] = QIter2(QQmwo,W6,Wo,dPj(va),wom(:,1),dPw(wom(:,3)),na); 
        [QQ.QQcwo] = QIter2(QQ.QQcwo,W6C,WoC,dPj(vc),WELL.WonC(:,1),dPw(WELL.WonC(:,3)),nc);
        [QQ.QQgwo] = QIter2(QQ.QQgwo,W6G,WoG,dPj(vg),WELL.WonG(:,1),dPw(WELL.WonG(:,3)),ng);
        [QQdwo] = QIter2(QQdwo,W6D,WoD,dPj(vd),wod(:,1),dPw(wod(:,3)),nd);    
        
        Fwater = AMW*Phj(1:na+nc+ng+nd,1) - dItimeW + [QQmwo(:,1) + qw;QQ.QQcwo(:,1);QQ.QQgwo(:,1);QQdwo(:,1) + qwd];
        
        Sw = Sw + (Fwater - SGMF.Cwp.*dPj )./SGMF.Cwsw;
        
        %Sw=Sw.*(Sw>=0).*(Sw<=1)+(Sw>1);
        So = 1 - Sw;
     
        Pw = Pw + dPw;
        Pj = PjOld + dPj;  
        Pa(v1==1)=Pj(va);
        Phia(v1==1,:)=Phj(va,:);
        Pd(v2==1)=Pj(vd);
        Phid(v2==1,:)=Phj(vd,:);
        Pc=Pj(vc);
        Pg=Pj(vg);
        Phic(:,:)=Phj(vc,:);
        Phig(:,:)=Phj(vg,:);
        flag_gim=sum(abs(dPj./Pj)>=1e-6)~=0;
        
    end
    
    [~,~,~,~,QQm,QQmwo]=Well_MKTF20(Pj(va),Pw(wom(:,3)),wom(:,2),wom(:,1),Uf(wom(:,3)),Cp(v_a(v1==1)),mu,CpW(wom(:,3)),kfw(va),kfo(va),SGMF,CMP);
    [~,~,~,~,QQ.QQc,QQ.QQcwo]=Well_MKTF20(Pj(vc),Pw(WELL.WonC(:,3)),WELL.WonC(:,2),WELL.WonC(:,1),Uf(WELL.WonC(:,3)),Cp(v_c),mu,CpW(WELL.WonC(:,3)),kfw(vc),kfo(vc),SGMF,CMP);
    [~,~,~,~,QQ.QQg,QQ.QQgwo]=Well_MKTF20(Pj(vg),Pw(WELL.WonG(:,3)),WELL.WonG(:,2),WELL.WonG(:,1),Uf(WELL.WonG(:,3)),Cp(v_g),mu,CpW(WELL.WonG(:,3)),kfw(vg),kfo(vg),SGMF,CMP);
    [~,~,~,~,QQd,QQdwo]=Well_MKTF20(Pj(vd),Pw(wod(:,3)),wod(:,2),wod(:,1),Uf(wod(:,3)),Cp(v_d(v2==1)),mu,CpW(wod(:,3)),kfw(vd),kfo(vd),SGMF,CMP);
    
    SW0([v_c,v_g])=Sw0([vc,vg]);
    SW([v_c,v_g])=Sw([vc,vg]);
    
    [BFRACM.Blc,BFRACM.Bwc,BFRACM.Boc]=crack_bond(Blc,Bwc,Boc,Pj(vc),Pj(va),A2CL,A2CW,A2CO,dt,ndt,c1);
    [BFRACM.Blg,BFRACM.Bwg,BFRACM.Bog]=crack_bond(Blg,Bwg,Bog,Pj(vg),Pj(va),A2GL,A2GW,A2GO,dt,ndt,c2);
 
    [BFRACM.Blcd,BFRACM.Bwcd,BFRACM.Bocd]=crack_bond(Blcd,Bwcd,Bocd,Pj(vc),Pj(vd),D2CL,D2CW,D2CO,dt,ndt,c1d);
    [BFRACM.Blgd,BFRACM.Bwgd,BFRACM.Bogd]=crack_bond(Blgd,Bwgd,Bogd,Pj(vg),Pj(vd),D2GL,D2GW,D2GO,dt,ndt,c2d);
    
    Pi=[Pa;Pc;Pg;Pd];
    Phi=[Phia;Phic;Phig;Phid];
    ndtI(i)=ndt;
    CL = sparse(1:nc,1:nc,CMPF.Cw(vc),nc,nc)*CW + CO;
    GL = sparse(1:ng,1:ng,CMPF.Cw(vg),ng,ng)*GW + GO;
  %  [ndt,j_ndt,fl]=vibor_t(ndt,Fl,Pc,Pg,Pw,PR,RC,dt,SW([v_c,v_g]),SW0([v_c,v_g]),CL,GL,dVc,dVg,W1C,W1G,...
  %      WELL.WonC,WELL.WonG,0,j_ndt);
    if fl==0
        fl2=fl2+1;
    end;
end;
%ndtI

%Bc(r1)=Bc;
%hj-1

%dSS=sum((Sw([vc,vg])-Sw0([vc,vg])).*[dVc;dVg])+sum(sQ1(:,1))+sum(Bwc);
%sum((SCwC(Na+1:end)-Sw0(Na+1:end)).*[dVc;dVg])

SW([v_a(v1==1),v_c,v_g,v_d(v2==1)])=Sw([va,vc,vg,vd]);

%CSw0=Sw0(vc);
%GSw=Sw(vg);
QQ.QQm(v1==1)=QQm; QQ.QQmwo(v1==1,:) = QQmwo;
QQ.QQd(v2==1)=QQd; QQ.QQdwo(v2==1,:) = QQdwo;

KWOG.w([v_a(v1==1),v_c,v_g,v_d(v2==1)]) = kfw;
KWOG.o([v_a(v1==1),v_c,v_g,v_d(v2==1)]) = kfo;

%KWOG.w([v_c,v_g])=kfw([vc,vg]);
%KWOG.o([v_c,v_g])=kfo([vc,vg]);
   
CMP.Bo([v_c,v_g],1) = CMPF.Bo([vc,vg],1);
CMP.Bw([v_c,v_g],1) = CMPF.Bw([vc,vg],1);
CMP.Mp([v_c,v_g],1) = CMPF.Mp([vc,vg],1);
   
end
function [v,vc,won,r1,c1,r2,c2,rc_gy,rc_in_h,T_gy,T_in_h]=ext_cr(CR_rc)
r1=CR_rc.r1;     c1=CR_rc.c1;
r2=CR_rc.r2;     c2=CR_rc.c2;

rc_gy=CR_rc.rc_gy;
rc_in_h=CR_rc.rc_in_h;

T_gy=CR_rc.T_gy;
T_in_h=CR_rc.T_in_h;
won=CR_rc.won;
v=CR_rc.v;
vc=CR_rc.vc;
end
function [Bl,Bw,Bo]=crack_bond(Bl,Bw,Bo,P1,P2,A2CL,A2CW,A2CO,dt,ndt,c)
 if isempty(c)==0
    Bw=Bw+(sum(A2CW,1)'.*P1-A2CW'*P2)*dt/ndt;
    Bo=Bo+(sum(A2CO,1)'.*P1-A2CO'*P2)*dt/ndt;
    Bl=Bl+(sum(A2CL,1)'.*P1-A2CL'*P2)*dt/ndt;
 end
end

