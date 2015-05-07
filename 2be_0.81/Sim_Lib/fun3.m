function [BFRACM,CSw,GSw,CSo,GSo,i,sQ1,sQ2,sQm,sQd,QWOG]=fun3(RC,Pi,Phi,SW,SO,Cp,PR,TRM,M2FR,WELL,Pw,dt,CR_rc,Qz,Qf,ndt,Pi0,GEOM,Lc,P0,SGM,ft,QWOG,VSAT)
                                                                              
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
MSw=SW(v_a);                                   DSw=SW(v_d); 
MSo=SO(v_a);                                   DSo=SO(v_d); 
MMp=SGM.Mp(v_a,:);                             DMp=SGM.Mp(v_d,:); 
MBwog=SGM.Bwog(v_a,:);                         DBwog=SGM.Bwog(v_d,:); 
MRs=SGM.Rs(v_a,:);                             DRs=SGM.Rs(v_d,:);

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

Pj(:,1)=[Pa(v1==1);Pc;Pg;Pd(v2==1);];
Phj=[Phia(v1==1,:);Phic;Phig;Phid(v2==1,:);];

Sw=[MSw(v1==1);SW(v_c);SW(v_g);DSw(v2==1)];
So=[MSo(v1==1);SO(v_c);SO(v_g);DSo(v2==1)];

SGMFR.Bwog=[MBwog(v1==1,:);SGM.Bwog(v_c,:);SGM.Bwog(v_g,:);DBwog(v2==1,:)];
SGMFR.Mp=[MMp(v1==1,:);SGM.Mp(v_c,:);SGM.Mp(v_g,:);DMp(v2==1,:)];
SGMFR.Rs=[MRs(v1==1,:);SGM.Rs(v_c,:);SGM.Rs(v_g,:);DRs(v2==1,:)];

Pgy=Pa(rc_gy(:,1));
Pgy2=Pd(rc_gy_d(:,1));
 
SO0=SO;
SW0=SW;

dPt = zeros(size(Sw));
dPA = zeros(size(Pa));
dPD = zeros(size(Pd));
dPw = zeros(size(Qz(:,1),1),1);

kfw=zeros(size(Sw));
kfo=zeros(size(Sw));
kfg=zeros(size(Sw));

KfwM=Sat_cal(MSw,1,1,as,aw); %water
KfoM=Sat_cal(MSw,2,1,as,aw); %oil
KfgM=Sat_cal(MSw,3,1,as,aw); %oil

KfwD=Sat_cal(DSw,1,1,ts,tw); %water
KfoD=Sat_cal(DSw,2,1,ts,tw); %oil
KfgD=Sat_cal(DSw,3,1,ts,tw); %oil

[bAl,bAw,bAo,ql,qw,qo]=Potok_GY0(T_gy,Pgy,Pa,rc_gy,KfwM,KfoM,KfgM,v1,mu,Na,vc1,SGM);
[bDl,bDw,bDo,qld,qwd,qod]=Potok_GY0(T_gy_d,Pgy2,Pd,rc_gy_d,KfwD,KfoD,KfgD,v2,mu,Nd,vc2,SGM);

fl=1;
i=0;
j_ndt=0;
fl2=0;
PjOld = Pj;
while fl2<2% & i<3000
    i = i + 1;
    kj = 0;
    flag_gim = 1;
    kfw(va)=Sat_cal(Sw(va),1,1,as,aw); %water
    kfo(va)=Sat_cal(Sw(va),2,1,as,aw); %oil
    kfg(va)=Sat_cal(Sw(va),3,1,as,aw);
    
    kfw([vc,vg,vd])=Sat_cal(Sw([vc,vg,vd]),1,1,ts,tw); %water
    kfo([vc,vg,vd])=Sat_cal(Sw([vc,vg,vd]),2,1,ts,tw); %oil
    kfg([vc,vg,vd])=Sat_cal(Sw([vc,vg,vd]),3,1,ts,tw);
    
    KfwM(v1==1)=kfw(va);      KfoM(v1==1)=kfo(va);      KfgM(v1==1)=kfg(va);
    KfwD(v2==1)=kfw(vd);      KfoD(v2==1)=kfo(vd);      KfgD(v2==1)=kfg(vd);
    
    KFA=[kfw(va),kfo(va),kfg(va)];
    KFC=[kfw(vc),kfo(vc),kfg(vc)];
    KFG=[kfw(vg),kfo(vg),kfg(vg)];
    KFD=[kfw(vd),kfo(vd),kfg(vd)];
    Sw0=Sw;
    So0=So;
    SGMFR.Rs(:,1) = SGMFR.Rs(:,2);
    SGMFR.Bwog(:,1) = SGMFR.Bwog(:,2);
    SGMFR.Bwog(:,3) = SGMFR.Bwog(:,4);
    SGMFR.Bwog(:,5) = SGMFR.Bwog(:,6);
    SGMFR.Mp(:,1) = SGMFR.Mp(:,2);
    
    while flag_gim==1 && kj<2
        kj=kj+1;
        [SGMFR]=SGim([dVa;dVc;dVg;dVd],Sw,So,SGMFR.Mp(:,2),zc,SGMFR.Rs(:,2),rs,SGMFR.Bwog,dPt,PjOld,1,dt/ndt,VSAT);
        PjOld = Pj;
        [vPa1,vPc1,vPg1,vPd1,dPa,dPc,dPg,dPd]=pre_potok_2(Phi,Phj,RC,rc_in_h,rc_in_hd,Na,Nd,na,nd);
        
        [~,TW,TO,TGG,TWCgs,TOCgs]=Potok_MKT_2(T_in,vPa1,Cp(v_a),mu,rc_in_h,Na,KfwM,KfoM,KfgM,kms(1),dPa,GEOM.L,Ro,Ke,SGMFR,v1);
        [CL,CW,CO,CG,CWCgs,COCgs]=Potok_Tube_2(TRM.TC,Pc,vPc1,kfw(vc),kfo(vc),kfg(vc),Cp(v_c),PR,RC.Cr2,RC.Cc2,kms(2),dPc,Lc,nc,SGMFR);
        [GL,GW,GO,GG,GWCgs,GOCgs]=Potok_Tube_2(TRM.TG,Pg,vPg1,kfw(vg),kfo(vg),kfg(vg),Cp(v_g),PR,RC.Gr2,RC.Gc2,kms(3),dPg,GEOM.Lg,ng,SGMFR);
        [~,DW,DO,DG,DWCgs,DOCgs]=Potok_MKT_2(TD_in,vPd1,Cp(v_d),mu,rc_in_hd,Nd,KfwD,KfoD,KfgD,kms(4),dPd,GEOM.L,Ro,Ke,SGMFR,v2);
        
        Pa1=Pi(v_a(v1==1));
        Cp1=Cp(v_a(v1==1));
        
        Pd1=Pi(v_d(v2==1));
        Cp1d=Cp(v_d(v2==1));
        
        [A2CL,A2CW,A2CO,A2CG,A2CWCgs,A2COCgs]=Obmen_T2M_2(TRM.TA2C,Pa1,Pc,mu,Cp1,Cp(v_c),r1,RC.ACc,KFA,KFC,SGMFR);
        [A2GL,A2GW,A2GO,A2GG,A2GWCgs,A2GOCgs]=Obmen_T2M_2(TRM.TA2G,Pa1,Pg,mu,Cp1,Cp(v_g),r2,RC.AGc,KFA,KFG,SGMFR);
        [A2DL,A2DW,A2DO,A2DG,A2DWCgs,A2DOCgs]=Obmen_T2M_2(TRM.TA2D,Pa1,Pd1,mu,Cp1,Cp1d,r3,c3,KFA,KFD,SGMFR);
        
        [D2CL,D2CW,D2CO,D2CG,D2CWCgs,D2COCgs]=Obmen_T2M_2(TRM.TD2C,Pd1,Pc,mu,Cp1d,Cp(v_c),r1d,RC.DCc,KFD,KFC,SGMFR);
        [D2GL,D2GW,D2GO,D2GG,D2GWCgs,D2GOCgs]=Obmen_T2M_2(TRM.TD2G,Pd1,Pg,mu,Cp1d,Cp(v_g),r2d,RC.DGc,KFD,KFG,SGMFR);
        
        [W1,W6,Wo,W7]=Well_MKT_2(wom(:,2),wom(:,1),Uf(wom(:,3)),Cp(v_a(v1==1)),mu,CpW(wom(:,3)),kfw(va),kfo(va),kfg(va),SGMFR);
        [W1C,W6C,WoC,W7C]=Well_MKT_2(WELL.WonC(:,2),WELL.WonC(:,1),Uf(WELL.WonC(:,3)),Cp(v_c),mu,CpW(WELL.WonC(:,3)),kfw(vc),kfo(vc),kfg(vc),SGMFR);
        [W1G,W6G,WoG,W7G]=Well_MKT_2(WELL.WonG(:,2),WELL.WonG(:,1),Uf(WELL.WonG(:,3)),Cp(v_g),mu,CpW(WELL.WonG(:,3)),kfw(vg),kfo(vg),kfg(vg),SGMFR);
        [W1D,W6D,WoD,W7D]=Well_MKT_2(wod(:,2),wod(:,1),Uf(wod(:,3)),Cp(v_d(v2==1)),mu,CpW(wod(:,3)),kfw(vd),kfo(vd),kfg(vd),SGMFR);
        
        AG1=TGG-sparse(1:na,1:na,sum(A2CG,2)+sum(A2GG,2)+sum(A2DG,2),na,na);
        CG1=CG-sparse(1:nc,1:nc,sum(A2CG,1)+sum(D2CG,1),nc,nc);
        GG1=GG-sparse(1:ng,1:ng,sum(A2GG,1)+sum(D2GG,1),ng,ng);
        DG1=DG-sparse(1:nd,1:nd,sum(A2DG,1)'+sum(D2CG,2)+sum(D2GG,2),nd,nd);
        
        AW1Cgs=TWCgs-sparse(1:na,1:na,sum(A2CWCgs,2)+sum(A2GWCgs,2)+sum(A2DWCgs,2),na,na);
        CW1Cgs=CWCgs-sparse(1:nc,1:nc,sum(A2CWCgs,1)+sum(D2CWCgs,1),nc,nc);
        GW1Cgs=GWCgs-sparse(1:ng,1:ng,sum(A2GWCgs,1)+sum(D2GWCgs,1),ng,ng);
        DW1Cgs=DWCgs-sparse(1:nd,1:nd,sum(A2DWCgs,1)'+sum(D2CWCgs,2)+sum(D2GWCgs,2),nd,nd);
        
        AO1Cgs=TOCgs-sparse(1:na,1:na,sum(A2COCgs,2)+sum(A2GOCgs,2)+sum(A2DOCgs,2),na,na);
        CO1Cgs=COCgs-sparse(1:nc,1:nc,sum(A2COCgs,1)+sum(D2COCgs,1),nc,nc);
        GO1Cgs=GOCgs-sparse(1:ng,1:ng,sum(A2GOCgs,1)+sum(D2GOCgs,1),ng,ng);
        DO1Cgs=DOCgs-sparse(1:nd,1:nd,sum(A2DOCgs,1)'+sum(D2COCgs,2)+sum(D2GOCgs,2),nd,nd);
        
        AW1=TW-sparse(1:na,1:na,sum(A2CW,2)+sum(A2GW,2)+sum(A2DW,2),na,na);
        CW1=CW-sparse(1:nc,1:nc,sum(A2CW,1)+sum(D2CW,1),nc,nc);
        GW1=GW-sparse(1:ng,1:ng,sum(A2GW,1)+sum(D2GW,1),ng,ng);
        DW1=DW-sparse(1:nd,1:nd,sum(A2DW,1)'+sum(D2CW,2)+sum(D2GW,2),nd,nd);
        
        AO1=TO-sparse(1:na,1:na,sum(A2CO,2)+sum(A2GO,2)+sum(A2DO,2),na,na);
        CO1=CO-sparse(1:nc,1:nc,sum(A2CO,1)+sum(D2CO,1),nc,nc);
        GO1=GO-sparse(1:ng,1:ng,sum(A2GO,1)+sum(D2GO,1),ng,ng);
        DO1=DO-sparse(1:nd,1:nd,sum(A2DO,1)'+sum(D2CO,2)+sum(D2GO,2),nd,nd);
        
        AMG=[AG1,   A2CG, A2GG, A2DG;
            A2CG',  CG1,  M2FR.C2G, D2CG';
            A2GG', M2FR.C2G', GG1,  D2GG';
            A2DG', D2CG, D2GL,  DG1];
        
        AMWCgs=[AW1Cgs,   A2CWCgs, A2GWCgs, A2DWCgs;
            A2CWCgs',  CW1Cgs,  M2FR.C2G, D2CWCgs';
            A2GWCgs', M2FR.C2G', GW1Cgs,  D2GWCgs';
            A2DWCgs', D2CWCgs, D2GWCgs,  DW1Cgs];
        
        AMOCgs=[AO1Cgs,   A2COCgs, A2GOCgs, A2DOCgs;
            A2COCgs',  CO1Cgs,  M2FR.C2G, D2COCgs';
            A2GOCgs', M2FR.C2G', GO1Cgs,  D2GOCgs';
            A2DOCgs', D2COCgs, D2GOCgs,  DO1Cgs];
        
        AMW=[AW1,   A2CW, A2GW, A2DW;
            A2CW',  CW1,  M2FR.C2G, D2CW';
            A2GW', M2FR.C2G', GW1,  D2GW';
            A2DW', D2CW, D2GW,  DW1];
        
        AMO=[AO1,   A2CO, A2GO, A2DO;
            A2CO',  CO1,  M2FR.C2G, D2CO';
            A2GO', M2FR.C2G', GO1,  D2GO';
            A2DO', D2CO, D2GO,  DO1];
        
        A1 = - sparse(1:na,1:na,SGMFR.Clp(va) + bAl',na,na)-sparse(wom(:,1),wom(:,1),W1,na,na);
        C1 = - sparse(1:nc,1:nc,SGMFR.Clp(vc)',nc,nc)-sparse(WELL.WonC(:,1),WELL.WonC(:,1),W1C,nc,nc);
        G1 = - sparse(1:ng,1:ng,SGMFR.Clp(vg)',ng,ng)-sparse(WELL.WonG(:,1),WELL.WonG(:,1),W1G,ng,ng);
        D1 = - sparse(1:nd,1:nd,SGMFR.Clp(vd) + bDl',nd,nd)-sparse(wod(:,1),wod(:,1),W1D,nd,nd);
        
        AM = [A1, zeros(na,nc), zeros(na,ng), zeros(na,nd);
            zeros(nc,na), C1, zeros(nc,ng), zeros(nc,nd);
            zeros(ng,na), zeros(ng,nc), G1, zeros(ng,nd);
            zeros(nd,na), zeros(nd,nc), zeros(nd,ng), D1];
        
        AM = AM + AMWCgs + AMOCgs + AMG;
        
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
        
        ba1 = ba1 - sparse(wom(:,1),ones(1,size(wom,1)),QWOG.QQm(wom(:,3)),na,1);
        bc1 = bc1 - sparse(WELL.WonC(:,1),ones(1,size(WELL.WonC,1)),QWOG.QQc,nc,1);
        bg1 = bg1 - sparse(WELL.WonG(:,1),ones(1,size(WELL.WonG,1)),QWOG.QQg,ng,1);
        bd1 = bd1 - sparse(wod(:,1),ones(1,size(wod,1)),QWOG.QQd(wod(:,3)),nd,1);
        
        [dItime,dItimeW,dItimeO] = NR_Time([dVa;dVc;dVg;dVd],Sw,So,Sw0,So0,SGMFR,dt/ndt);
        
        dPhi = AMWCgs*Phj(:,1) + AMOCgs*Phj(:,2) + AMG*Phj(:,3);
        
        %   BC = [ba1+ql-bAl'.*dPt(rc_gy(:,1));bc1;bg1;bd1+qld-bDl'.*dPt(rc_gy_d(:,1))] - dItime + dPhi;
        
        BC = [ba1 - ql;bc1;bg1;bd1 - qld] + dItime - dPhi;
        
        Qz1 = sparse(wom(Qf~=0,3),ones(1,size(wom(Qf~=0,3),1)),Qz(wom(Qf~=0,3),ft+1),nw,1);
        Qm1 = sparse(wom(Qf~=0,3),ones(1,size(wom(Qf~=0,3),1)),QWOG.QQm(wom(Qf~=0,3),ft+1),nw,1);
        Qc1 = sparse(WELL.WonC(Qf1(WELL.WonC(:,3))~=0,3),ones(1,size(WELL.WonC(Qf1(WELL.WonC(:,3))~=0,3),1)),QWOG.QQc(Qf1(WELL.WonC(:,3))~=0),nw,1);
        Qg1 = sparse(WELL.WNG(Qf1(WELL.WonG(:,3))~=0),ones(1,size(WELL.WonG(Qf1(WELL.WonG(:,3))~=0,3),1)),QWOG.QQg(Qf1(WELL.WonG(:,3))~=0),nw,1);
        Qd1 = sparse(wod(Qf(wod(:,3))~=0,3),ones(1,size(wod(Qf(wod(:,3))~=0,3),1)),QWOG.QQd(wod(Qf(wod(:,3)),3)~=0),nw,1);
        
        Q11 = -Qz1+Qm1+Qc1+Qg1+Qd1;
        
        dPt=[BC',Q11(wom(:,3))']/[AM,WM2;WM1,WM3];
        dPt = dPt';
        dPj(:,1)=dPt(1:na+nc+ng+nd);
        dPA(v1==1)=dPj(va);
        dPD(v2==1)=dPj(vd);
        dPw(wom(Qf~=0,3))=dPt(na+nc+ng+nd+1:end);
        
        [bAl,bAw,bAo,ql,qw,qo]=Potok_GY(T_gy,Pgy,Pa,dPA(v_a),rc_gy,KfwM,KfoM,KfgM,v1,mu,Na,vc1,SGMFR,ql,qw,qo);
        [bDl,bDw,bDo,qld,qwd,qod]=Potok_GY(T_gy_d,Pgy2,Pd,dPD(v_d),rc_gy_d,KfwD,KfoD,KfgD,v2,mu,Nd,vc2,SGMFR,qld,qwd,qod);
        
        [QWOG.QQm,QQmwo]=QIter(QWOG.QQm,W1,W6,Wo,dPj(va),wom(:,1),dPw(wom(:,3)),SGMFR);
        [QWOG.QQc,QQcwo]=QIter(QWOG.QQc,W1C,W6C,WoC,dPj(vc),WELL.WonC(:,1),dPw(WELL.WonC(:,3)),SGMFR);
        [QWOG.QQg,QQgwo]=QIter(QWOG.QQg,W1G,W6G,WoG,dPj(vg),WELL.WonG(:,1),dPw(WELL.WonG(:,3)),SGMFR);
        [QWOG.QQd,QQdwo]=QIter(QWOG.QQd,W1D,W6D,WoD,dPj(vd),wod(:,1),dPw(wod(:,3)),SGMFR);
        
        QWOG.QQmwo(v1==1,1) = QWOG.QQmwo(v1==1,1) + sparse(wom(:,1),ones(1,size(wom,1)),QQmwo(:,1),na,1);
        QWOG.QQcwo(:,1) = QWOG.QQcwo(:,1) + sparse(WELL.WonC(:,1),ones(1,size(WELL.WonC,1)),QQcwo(:,1),nc,1);
        QWOG.QQgwo(:,1) = QWOG.QQgwo(:,1) + sparse(WELL.WonG(:,1),ones(1,size(WELL.WonG,1)),QQgwo(:,1),ng,1);
        QWOG.QQdwo(v2==1,1) = QWOG.QQdwo(v2==1,1) + sparse(wod(:,1),ones(1,size(WELL.WonD,1)),QQdwo(:,1),nd,1);
        
        QWOG.QQmwo(v1==1,2) = QWOG.QQmwo(v1==1,2) + sparse(wom(:,1),ones(1,size(wom,1)),QQmwo(:,2),na,1);
        QWOG.QQcwo(:,2) = QWOG.QQcwo(:,2) + sparse(WELL.WonC(:,1),ones(1,size(WELL.WonC,1)),QQcwo(:,2),nc,1);
        QWOG.QQgwo(:,2) = QWOG.QQgwo(:,2) + sparse(WELL.WonG(:,1),ones(1,size(WELL.WonG,1)),QQgwo(:,2),ng,1);
        QWOG.QQdwo(v2==1,2) = QWOG.QQdwo(v2==1,2) + sparse(wod(:,1),ones(1,size(WELL.WonD,1)),QQdwo(:,2),nd,1);
        
        Phj(1:na+nc+ng+nd,1) = Phj(1:na+nc+ng+nd,1) + dPj;   %//////////////
        Phj(1:na+nc+ng+nd,2) = Phj(1:na+nc+ng+nd,2) + dPj;
        
        Fwater = AMW*Phj(1:na+nc+ng+nd,1) - dItimeW + [QWOG.QQmwo(v1==1,1) + qw;QWOG.QQcwo(:,1);QWOG.QQgwo(:,1);QWOG.QQdwo(v2==1,1) + qwd];
        Foil = AMO*Phj(1:na+nc+ng+nd,2) - dItimeO + [QWOG.QQmwo(v1==1,2) + qo;QWOG.QQcwo(:,2);QWOG.QQgwo(:,2);QWOG.QQdwo(v2==1,1) + qod];
        
        Sw = Sw + (Fwater - SGMFR.Cwp.*dPt(1:na+nc+ng+nd))./SGMFR.Cws;
        Sw=Sw.*(Sw>=0).*(Sw<=1)+(Sw>1);
        So = 1 - Sw;
        %So = So + (Foil - SGMFR.Cop*dPt(1:na+nc+ng+nd))./SGMFR.Cos;
        
        Pj = PjOld + dPj;     %////////////////
        Pa(v1==1)=Pj(va);
        Pd(v2==1)=Pj(vd);
        Pc=Pj(vc);
        Pg=Pj(vg);
        
        flag_gim=sum(abs(dPj./Pj)>=1e-6)~=0;
    end
    
    SW0([v_c,v_g])=Sw0([vc,vg]);
    SW([v_c,v_g])=Sw([vc,vg]);
    
    [BFRACM.Blc,BFRACM.Bwc,BFRACM.Boc]=crack_bond(Blc,Bwc,Boc,Pj(vc),Pj(va),A2CL,A2CW,A2CO,dt,ndt,c1);
    [BFRACM.Blg,BFRACM.Bwg,BFRACM.Bog]=crack_bond(Blg,Bwg,Bog,Pj(vg),Pj(va),A2GL,A2GW,A2GO,dt,ndt,c2);
    
    [BFRACM.Blcd,BFRACM.Bwcd,BFRACM.Bocd]=crack_bond(Blcd,Bwcd,Bocd,Pj(vc),Pj(vd),D2CL,D2CW,D2CO,dt,ndt,c1d);
    [BFRACM.Blgd,BFRACM.Bwgd,BFRACM.Bogd]=crack_bond(Blgd,Bwgd,Bogd,Pj(vg),Pj(vd),D2GL,D2GW,D2GO,dt,ndt,c2d);
    
    Pi=[Pa;Pc;Pg;Pd];
    ndtI(i)=ndt;
    [ndt,j_ndt,fl]=vibor_t(ndt,Fl,Pc,Pg,Pw,PR,RC,dt,SW([v_c,v_g]),SW0([v_c,v_g]),CL,GL,dVc,dVg,W1C,W1G,...
        WonC,WonG,0,j_ndt);
    if fl==0
        fl2=fl2+1;
    end;
    
    sQm(:,:)=sQm + QBild(QWOG.QQm(wom(:,3)),QWOG.QQmwo(WELL.Won(wom(:,3)),:),Uf(wom(:,3)),wom(:,1),dt/ndt,wom(:,3),nw,W1);
    sQ1(:,:)=sQ1 + QBild(QWOG.QQc,QWOG.QQcwo(WELL.WonC(:,1),:),Uf(WELL.WonC(:,3)),WELL.WonC(:,1),dt/ndt,WELL.WonC(:,3),nw,W1C);
    sQ2(:,:)=sQ2 + QBild(QWOG.QQg,QWOG.QQgwo(WELL.WonG(:,1),:),Uf(WELL.WonG(:,3)),WELL.WonG(:,1),dt/ndt,WELL.WonG(:,3),nw,W1G);
    sQd(:,:)=sQd + QBild(QWOG.QQd(wod(:,3)),QWOG.QQdwo(WELL.Won(wod(:,3)),:),Uf(wod(:,3)),wod(:,1),dt/ndt,wod(:,3),nw,W1D);
end;
%ndtI

%Bc(r1)=Bc;
%hj-1

%dSS=sum((Sw([vc,vg])-Sw0([vc,vg])).*[dVc;dVg])+sum(sQ1(:,1))+sum(Bwc);
%sum((SCwC(Na+1:end)-Sw0(Na+1:end)).*[dVc;dVg])

CSw=Sw(vc); CSo=So(vc);
GSw=Sw(vg); GSo=So(vg); 
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
    Bw=Bw+(A2CW*P1-sum(A2CW,2).*P2)*dt/ndt;
    Bo=Bo+(A2CO*P1-sum(A2CO,2).*P2)*dt/ndt;
    Bl=Bl+(A2CL*P1-sum(A2CL,2).*P2)*dt/ndt;
 end
end

