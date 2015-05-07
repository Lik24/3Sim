function [Bwc,Bwg,Boc,Bog,Bwcd,Bwgd,Bocd,Bogd,CSw,GSw,i,Q1,Q2,Qm,Qd,dSS,ndt]=fun1(RC,Pi,SW,Cp,PR,TC,TG,a2c,a2g,a2d,d2c,d2g,WonC,...
    WonG,Uf,CpW,Pw,dt,dV,CR_rc,Qz,Qf,ndt,Pi0,L,Lc,Lg,Ke,dZA,Mp,Bwo,P0)

as=PR.as;
aw=PR.aw;
ts=PR.ts;
tw=PR.tw;
mu=PR.mu;
kms=PR.kms;
Ro=PR.Ro;
zc=PR.zc;

Na=RC.na;
nc=RC.nc;
ng=RC.ng;
Nd=RC.nd;
C2G=sparse(nc,ng);     C2GL=C2G;  C2GW=C2G;

[v1,vc1,wom,r1,c1,r2,c2,rc_gy,rc_in_h,T_gy,T_in]=ext_cr(CR_rc(1,1));
[v2,vc2,wod,r1d,c1d,r2d,c2d,rc_gy_d,rc_in_hd,T_gy_d,TD_in]=ext_cr(CR_rc(1,2));
a2d=CR_rc(1,3).a2d;
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

Pa=Pi(v_a);    Pc=Pi(v_c);    Pg=Pi(v_g);      Pd=Pi(v_d);
dVa=dV(v_a);   dVc=dV(v_c);   dVg=dV(v_g);     dVd=dV(v_d);
MSw=SW(v_a);                                   DSw=SW(v_d);
MMp=Mp(v_a,:);                                 DMp=Mp(v_d,:);
MBwo=Bwo(v_a,:);                               DBwo=Bwo(v_d,:);

dVa(v1~=1)=[];   dVd(v2~=1)=[];

Qz=Qz/dt;
nw=size(Qz,1);
Qf=Qf(unique(wom(:,3)));

%%
vP1=Pa(RC.ACr)>=Pc(RC.ACc);
vP2=Pa(RC.AGr)>=Pg(RC.AGc);
fl1=sum([SW(v_c);SW(v_g)])/size([SW(v_c);SW(v_g)],1)==1;
fl2=sum([vP1;vP2])==0;
Fl=fl1*fl2;

%%

Bwc=zeros(size(RC.ACr,1),1);   Boc=Bwc;
Bwg=zeros(size(RC.AGr,1),1);   Bog=Bwg;

Bwcd=zeros(size(RC.DCr,1),1);   Bocd=Bwcd;
Bwgd=zeros(size(RC.DGr,1),1);   Bogd=Bwgd;

Q1=zeros(nw,5);
Q2=zeros(nw,5);
Qm=zeros(nw,5);
Qd=zeros(nw,5);

Pj(:,1)=[Pa(v1==1);Pc;Pg;Pd(v2==1);];
Sw=[MSw(v1==1);SW(v_c);SW(v_g);DSw(v2==1)];
Mp=[MMp(v1==1,:);Mp(v_c,:);Mp(v_g,:);DMp(v2==1,:)];
Bwo=[MBwo(v1==1,:);Bwo(v_c,:);Bwo(v_g,:);DBwo(v2==1,:)];

Pgy=Pa(rc_gy(:,1));
Pgy2=Pd(rc_gy_d(:,1));

Sw0=Sw;
SW0=SW;

kfw=zeros(size(Sw));
kfo=zeros(size(Sw));

KfwM=Sat_cal(MSw,1,1,as,aw); %water
KfoM=Sat_cal(MSw,2,1,as,aw); %oil

KfwD=Sat_cal(DSw,1,1,ts,tw); %water
KfoD=Sat_cal(DSw,2,1,ts,tw); %oil

wna=unique(wom(:,3));
wnc=unique(WonC(:,3));
wng=unique(WonG(:,3));
wnd=unique(wod(:,3));

i=0;
j_ndt=1/ndt;
fl2=0;
DZA=pre2Pot(dZA,[na,nc,ng,nd],RC,rc_in_h,rc_in_hd);
bicflag=1;
Pj0=Pi0([v_a(v1==1),v_c,v_g,v_d(v2==1)]);
while fl2<2% & i<3000
    i=i+1;
    
    kfw0(va,1)=Sat_cal(Sw(va),1,1,as,aw); %water
    kfo0(va,1)=Sat_cal(Sw(va),2,1,as,aw); %oil
    
    kfw0([vc,vg,vd],1)=Sat_cal(Sw([vc,vg,vd]),1,1,ts,tw); %water
    kfo0([vc,vg,vd],1)=Sat_cal(Sw([vc,vg,vd]),2,1,ts,tw); %oil
    
    [vPa1,vPc1,vPg1,vPd1,dPa,dPc,dPg,dPd]=pre_potok_2(Pi,Pj,RC,rc_in_h,rc_in_hd,Na,Nd,na,nd,dZA,DZA);
    
    kj=0;
    flag_gim=1;
  
     Pt0=Pj;
    while flag_gim==1 && kj<20
        kj=kj+1;
        
        [Clp,Cwp,Cws,A,Bwo,Mp]=SGim([dVa;dVc;dVg;dVd],Sw,Mp,zc,Bwo,Pi,1,P0([v_a(v1==1),v_c,v_g,v_d(v2==1)]),va,vc,vg,vd,zeros(1,0),dt/ndt);
        
        kfw=kfw0./Bwo(:,2);           kfo=kfo0./Bwo(:,4);
        
        KfwM(v1==1)=kfw(va);      KfoM(v1==1)=kfo(va);
        KfwD(v2==1)=kfw(vd);      KfoD(v2==1)=kfo(vd);
        
        [TL2,TW2]=Potok_MKT_2(T_in,vPa1,Cp(v_a),mu,rc_in_h,Na,KfwM,KfoM,kms(1),dPa,L,Ro,Ke,A(va),dZ(1,:),v1);
        [CL1,CW1,CL,CW]=Potok_Tube_2(TC,Pc,vPc1,kfw(vc),kfo(vc),Cp(v_c),PR,RC.Cr2,RC.Cc2,kms(2),dPc,Lc,nc,A(vc),dZ(2,:));
        [GL1,GW1,GL,GW]=Potok_Tube_2(TG,Pg,vPg1,kfw(vg),kfo(vg),Cp(v_g),PR,RC.Gr2,RC.Gc2,kms(3),dPg,Lg,ng,A(vg),dZ(3,:));
        [DL2,DW2]=Potok_MKT_2(TD_in,vPd1,Cp(v_d),mu,rc_in_hd,Nd,KfwD,KfoD,kms(4),dPd,L,Ro,Ke,A(vd),dZ(4,:),v2);
        
        Pa1=Pi(v_a(v1==1));
        Cp1=Cp(v_a(v1==1));
        
        Pd1=Pi(v_d(v2==1));
        Cp1d=Cp(v_d(v2==1));
        
        KFA=[kfw(va),kfo(va)];
        KFC=[kfw(vc),kfo(vc)];
        KFG=[kfw(vg),kfo(vg)];
        KFD=[kfw(vd),kfo(vd)];

        A2C=Obmen_T2M_2(a2c,Pa1,Pc,mu,Cp1,Cp(v_c),r1,RC.ACc,KFA,KFC,A(va),A(vc));
        A2G=Obmen_T2M_2(a2g,Pa1,Pg,mu,Cp1,Cp(v_g),r2,RC.AGc,KFA,KFG,A(va),A(vg));
        A2D=Obmen_T2M_2(a2d,Pa1,Pd1,mu,Cp1,Cp1d,r3,c3,KFA,KFD,A(va),A(vd));
        
        D2C=Obmen_T2M_2(d2c,Pd1,Pc,mu,Cp1d,Cp(v_c),r1d,RC.DCc,KFD,KFC,A(vd),A(vc));
        D2G=Obmen_T2M_2(d2g,Pd1,Pg,mu,Cp1d,Cp(v_g),r2d,RC.DGc,KFD,KFG,A(vd),A(vg));
        
        [W1,W6,W7]=Well_MKT_2(wom(:,2),wom(:,1),Uf(wom(:,3)),Cp(v_a(v1==1)),mu,CpW(wom(:,3)),kfw(va),kfo(va),A(va));
        [W1C,W6C,W7C]=Well_MKT_2(WonC(:,2),WonC(:,1),Uf(WonC(:,3)),Cp(v_c),mu,CpW(WonC(:,3)),kfw(vc),kfo(vc),A(vc));
        [W1G,W6G,W7G]=Well_MKT_2(WonG(:,2),WonG(:,1),Uf(WonG(:,3)),Cp(v_g),mu,CpW(WonG(:,3)),kfw(vg),kfo(vg),A(vg));
        [W1D,W6D,W7D]=Well_MKT_2(wod(:,2),wod(:,1),Uf(wod(:,3)),Cp(v_d(v2==1)),mu,CpW(wod(:,3)),kfw(vd),kfo(vd),A(vd));
        
        
        C1=CL1-sparse(1:nc,1:nc,sum(A2C.L2,1)+sum(D2C.L2,1)+Clp(vc)',nc,nc)-sparse(WonC(:,1),WonC(:,1),W1C,nc,nc);
        G1=GL1-sparse(1:ng,1:ng,sum(A2G.L2,1)+sum(D2G.L2,1)+Clp(vg)',ng,ng)-sparse(WonG(:,1),WonG(:,1),W1G,ng,ng);
        
        C2=CW1-sparse(1:nc,1:nc,sum(A2C.W,1)+sum(D2C.W,1)+Cwp(vc)',nc,nc)-sparse(WonC(:,1),WonC(:,1),W6C,nc,nc);
        G2=GW1-sparse(1:ng,1:ng,sum(A2G.W,1)+sum(D2G.W,1)+Cwp(vg)',ng,ng)-sparse(WonG(:,1),WonG(:,1),W6G,ng,ng);
        
        [bAl,bAw,bl,bw]=Potok_GY(T_gy,Pgy,Pi0(v_a),rc_gy,KfwM,KfoM,vc1,mu,Na,A(va));
        [bDl,bDw,bld,bwd]=Potok_GY(T_gy_d,Pgy2,Pd,rc_gy_d,KfwD,KfoD,vc2,mu,Nd,A(vd));
        
        [Gr,Grw]=Gravity(TL2,TW2,CL1,CW1,GL1,GW1,DL2,DW2,[],A,dZ);
        
        A1=TL2-sparse(1:na,1:na,sum(A2C.L1,2)+sum(A2G.L1,2)+sum(A2D.L1,2)+Clp(va)+bAl',na,na)-sparse(wom(:,1),wom(:,1),W1,na,na);
        D1=DL2-sparse(1:nd,1:nd,sum(D2C.L1,2)+sum(D2G.L1,2)+sum(A2D.L2,1)'+Clp(vd)+bDl',nd,nd)-sparse(wod(:,1),wod(:,1),W1D,nd,nd);
        
        AMC1=[A1,   A2C.L1,  A2G.L1,  A2D.L1;
                   A2C.L2',  C1,   C2GL,  D2C.L2';
                   A2G.L2', C2GL',  G1,   D2G.L2';
                   A2D.L2', D2C.L1,  D2G.L1,   D1];
        
        ba1=sparse(wom(:,1),ones(1,size(wom,1)),-W1.*Pw(wom(:,3)),na,1);
        bc1=sparse(WonC(:,1),ones(1,size(WonC,1)),-W1C.*Pw(WonC(:,3)),nc,1);
        bg1=sparse(WonG(:,1),ones(1,size(WonG,1)),-W1G.*Pw(WonG(:,3)),ng,1);
        bd1=sparse(wod(:,1),ones(1,size(wod,1)),-W1D.*Pw(wod(:,3)),nd,1);
        
        W2M=sparse(wom(:,3),wom(:,1),W1,nw,na);
        W2C=sparse(WonC(:,3),WonC(:,1),W1C,nw,nc);
        W2G=sparse(WonG(:,3),WonG(:,1),W1G,nw,ng);
        W2D=sparse(wod(:,3),wod(:,1),W1D,nw,nd);
        
        ba1=ba1.*(sum(W2M(Qz==0,:),1)~=0)';
        bc1=bc1.*(sum(W2C(Qz==0,:),1)~=0)';
        bg1=bg1.*(sum(W2G(Qz==0,:),1)~=0)';
        bd1=bd1.*(sum(W2D(Qz==0,:),1)~=0)';
        
        WM1=[W2M,W2C,W2G,W2D];
        WM2=WM1';
        W3vec=sparse(wom(:,3),1,W1,nw,1)+sparse(WonC(:,3),1,W1C,nw,1)+sparse(WonG(:,3),1,W1G,nw,1)+sparse(wod(:,3),1,W1D,nw,1);
        WM3=-sparse(1:nw,1:nw,W3vec,nw,nw);
        %WM1(Qz~=0,:)
        % full([Qz(wn(Qf~=0)),W3vec(wn(Qf~=0)),Pw(wn(Qf~=0))])
        
        WM1=WM1(wna(Qf~=0),:);
        WM2=WM2(:,wna(Qf~=0));
        WM3=WM3(wna(Qf~=0),wna(Qf~=0));
        
        BC1=[ba1-bl';bc1;bg1;bd1-bld']'-(Clp.*Pj0)'-Gr';
        
        AA=[AMC1,WM2;WM1,WM3];
        BB=[BC1,Qz(wna(Qf~=0))']';
%         if bicflag==1
%             [L,U] = ilu(AA,struct('type','ilutp','droptol',1e-5));
%         end
%         [Pt,~] = bicgstab(AA,BB,10e-5,2000,L,U,[Pj;Pw(wna(Qf~=0))]);
        Pt(:,1)=BB'/AA;
        %Pt=[BC1,Qz(wna(Qf~=0))']/[AMC1,WM2;WM1,WM3];
       
        Pj(:,1)=Pt(1:na+nc+ng+nd);
        
        flag_gim=sum(abs(Pj-Pt0(1:na+nc+ng+nd))./Pj>=1e-6)~=0;
        flag_gim=flag_gim*(sum(zc==0)==0);
       % flag_gim=0;
       Pt0=Pt;
    end
    
    Pj0=Pj;
    Pw(wna(Qf~=0))=Pt(na+nc+ng+nd+1:end);
    
    A2=TW2-sparse(1:na,1:na,sum(A2C.W,2)+sum(A2G.W,2)+sum(A2D.W,2)+Cwp(va)+bAw',na,na)-sparse(wom(:,1),wom(:,1),W6,na,na);
    D2=DW2-sparse(1:nd,1:nd,sum(D2C.W,2)+sum(D2G.W,2)+sum(A2D.W,1)'+Cwp(vd)+bDw',nd,nd)-sparse(wod(:,1),wod(:,1),W6D,nd,nd);
    
    AMC2=[A2,   A2C.W,  A2G.W,  A2D.W;
            A2C.W',  C2,   C2GW,  D2C.W';
            A2G.W', C2GW',  G2,   D2G.W';
            A2D.W', D2C.W,  D2G.W,   D2];
    
    ba2=sparse(wom(:,1),ones(1,size(wom,1)),W6.*Pw(wom(:,3)),na,1);
    bc2=sparse(WonC(:,1),ones(1,size(WonC,1)),W6C.*Pw(WonC(:,3)),nc,1);
    bg2=sparse(WonG(:,1),ones(1,size(WonG,1)),W6G.*Pw(WonG(:,3)),ng,1);
    bd2=sparse(wod(:,1),ones(1,size(wod,1)),W6D.*Pw(wod(:,3)),nd,1);
    
    ba2=ba2+bw';
    bd2=bd2+bwd';
    BC=[ba2;bc2;bg2;bd2]+(Cwp.*Pj)+Grw;
    
    Pc=Pj(vc);
    Pg=Pj(vg);
    
    Sw1=Sw;
    Sw=Sw1+(AMC2*Pj+BC)./Cws*dt/ndt;
    Sw=Sw.*(Sw>=0).*(Sw<=1)+(Sw>1);
    %fghgfh
    bicflag=sum(abs(Sw-Sw1)>0.05)>1;
    
    SW0([v_c,v_g])=Sw1([vc,vg]);
    SW([v_c,v_g])=Sw([vc,vg]);
    %hj(:,i)=SCw(vc);
    
    [Boc,Bwc]=crack_bond(Boc,Bwc,Pj(vc),Pj(va),A2C,dt,ndt,c1);
    [Bog,Bwg]=crack_bond(Bog,Bwg,Pj(vg),Pj(va),A2G,dt,ndt,c2);
    
    [Bocd,Bwcd]=crack_bond(Bocd,Bwcd,Pj(vc),Pj(vd),D2C,dt,ndt,c1d);
    [Bogd,Bwgd]=crack_bond(Bogd,Bwgd,Pj(vg),Pj(vd),D2G,dt,ndt,c2d);
    
    % full([Pj(won)<Pw(wn),Uf(wn)])
    % dfgh
    
    Pa(v1==1)=Pj(va);
    Pd(v2==1)=Pj(vd);
    Pi=[Pa;Pc;Pg;Pd];
    
    Qm(:,:)=Qm+QBild2(wom(:,2),wom(:,1),kfw(va).*Bwo(va,2),kfo(va).*Bwo(va,4),PR.mu,Uf(wom(:,3)),Pj(va,1),Pw(wom(:,3)),dt/ndt,wom(:,3),nw,wna);
    Q1(:,:)=Q1+QBild2(WonC(:,2),WonC(:,1),kfw(vc).*Bwo(vc,2),kfo(vc).*Bwo(vc,4),PR.mu,Uf(WonC(:,3)),Pj(vc,1),Pw(WonC(:,3)),dt/ndt,WonC(:,3),nw,wnc);
    Q2(:,:)=Q2+QBild2(WonG(:,2),WonG(:,1),kfw(vg).*Bwo(vg,2),kfo(vg).*Bwo(vg,4),PR.mu,Uf(WonG(:,3)),Pj(vg,1),Pw(WonG(:,3)),dt/ndt,WonG(:,3),nw,wng);
    Qd(:,:)=Qd+QBild2(wod(:,2),wod(:,1),kfw(vd).*Bwo(vd,2),kfo(vd).*Bwo(vd,4),PR.mu,Uf(wod(:,3)),Pj(vd,1),Pw(wod(:,3)),dt/ndt,wod(:,3),nw,wnd);
    
    %SCwC(Na+1:Na+nc)-1
    ndtI(i)=ndt;
    [ndt,j_ndt,fl]=vibor_t(ndt,Fl,Pc,Pg,Pw,PR,RC,dt,Sw([vc,vg]),Sw1([vc,vg]),CL,GL,dVc,dVg,W1C,W1G,...
        WonC,WonG,0,j_ndt);
    if fl==0
        fl2=fl2+1;
        %     if fl2==2
        %     ndt=ndtI(end-1);
        %     end
    end;
    % [Pj(va),Pj(vc(end:-1:1)),Pj(vd)]
end;
ndt=ndtI(end);

%Bc(r1)=Bc;
%hj-1

dSS=sum((Sw([vc,vg])-Sw0([vc,vg])).*[dVc;dVg])+sum(Q1(:,1))+sum(Bwc);
%sum((SCwC(Na+1:end)-Sw0(Na+1:end)).*[dVc;dVg])

CSw=Sw(vc);
GSw=Sw(vg);
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
function [Bo,Bw]=crack_bond(Bo,Bw,P1,P2,A2C,dt,ndt,c)
if isempty(c)==0
    Bw=Bw+(sum(A2C.W,1)'.*P1-A2C.W'*P2)*dt/ndt;
    Bo=Bo+(sum(A2C.O,1)'.*P1-A2C.O'*P2)*dt/ndt;
end
end

