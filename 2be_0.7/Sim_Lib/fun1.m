function [Bwc,Bwg,Blc,Blg,Bwcd,Bwgd,Blcd,Blgd,CSw,GSw,i,Q1,Q2,Qm,Qd,dSS]=fun1(RC,Pi,SW,Cp,PR,TC,TG,A2C,A2G,A2D,D2C,D2G,WonC,...
    WonG,Uf,CpW,Pw,dt,dV,CR_rc,Qz,Qf,ndt,Pi0,L,Lc,Lg,Ke)

as=PR.as;
aw=PR.aw;
ts=PR.ts;
tw=PR.tw;
mu=PR.mu;
kms=PR.kms;
Ro=PR.Ro;

Na=RC.na;
nc=RC.nc;
ng=RC.ng;
Nd=RC.nd;
C2G=sparse(nc,ng);     C2GL=C2G;  C2GW=C2G;

[v1,wom,r1,c1,r2,c2,rc_gy,rc_in_h,T_gy,T_in]=ext_cr(CR_rc(1,1));
[v2,wod,r1d,c1d,r2d,c2d,rc_gy_d,rc_in_hd,T_gy_d,TD_in]=ext_cr(CR_rc(1,2));
A2D=CR_rc(1,3).a2d;
r3=CR_rc(1,3).r;    c3=CR_rc(1,3).c;

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
dVa(v1~=1)=[];   dVd(v2~=1)=[];
 
Qz=Qz/dt;
nw=size(Qz,1);
Qf=Qf(wom(:,3));

%%
vP1=Pa(RC.ACr)>=Pc(RC.ACc);
vP2=Pa(RC.AGr)>=Pg(RC.AGc);
    fl1=sum([SW(v_c);SW(v_g)])/size([SW(v_c);SW(v_g)],1)==1;
    fl2=sum([vP1;vP2])==0;
    Fl=fl1*fl2;
    
%%

Bwc=zeros(size(A2C,1),1);   Blc=Bwc;
Bwg=zeros(size(A2G,1),1);   Blg=Bwg;

Bwcd=zeros(size(D2C,1),1);   Blcd=Bwcd;
Bwgd=zeros(size(D2G,1),1);   Blgd=Bwgd;

Q1=zeros(nw,5);
Q2=zeros(nw,5);
Qm=zeros(nw,5);
Qd=zeros(nw,5);

Pj(:,1)=[Pa(v1==1);Pc;Pg;Pd(v2==1);];
Sw=[MSw(v1==1);SW(v_c);SW(v_g);DSw(v2==1)];

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

fl=1;
i=0;
j_ndt=0;

while fl==1% & i<3000
i=i+1;

     kfw(va)=Sat_cal(Sw(va),1,1,as,aw); %water
     kfo(va)=Sat_cal(Sw(va),2,1,as,aw); %oil
     
     kfw([vc,vg,vd])=Sat_cal(Sw([vc,vg,vd]),1,1,ts,tw); %water
     kfo([vc,vg,vd])=Sat_cal(Sw([vc,vg,vd]),2,1,ts,tw); %oil
     
     KfwM(v1==1)=kfw(va);      KfoM(v1==1)=kfo(va);
     KfwD(v2==1)=kfw(vd);      KfoD(v2==1)=kfo(vd);
     
     [vPa1,vPc1,vPc2,vPg1,vPg2,vPd1,dPa,dPc,dPg,dPd]=pre_potok_2(Pi,Pj,RC,rc_in_h,rc_in_hd,Na,Nd,na,nd);
     
     [TL2,TW2]=Potok_MKT_2(T_in,vPa1,Cp(v_a),mu,rc_in_h,Na,KfwM,KfoM,kms(1),dPa,L,Ro,Ke);
     [CL,CW]=Potok_Tube_2(TC,Pc,vPc1,vPc2,kfw(vc),kfo(vc),Cp(v_c),PR,RC.Cr2,RC.Cc2,kms(2),dPc,Lc,nc);
     [GL,GW]=Potok_Tube_2(TG,Pg,vPg1,vPg2,kfw(vg),kfo(vg),Cp(v_g),PR,RC.Gr2,RC.Gc2,kms(3),dPg,Lg,ng);
     [DL2,DW2]=Potok_MKT_2(TD_in,vPd1,Cp(v_d),mu,rc_in_hd,Nd,KfwD,KfoD,kms(4),dPd,L,Ro,Ke);

     Pa1=Pi(v_a(v1==1));
     Cp1=Cp(v_a(v1==1));
     
     Pd1=Pi(v_d(v2==1));
     Cp1d=Cp(v_d(v2==1));

     KFA=[kfw(va),kfo(va)];
     KFC=[kfw(vc),kfo(vc)];
     KFG=[kfw(vg),kfo(vg)];
     KFD=[kfw(vd),kfo(vd)];
     
     [A2CL,A2CW]=Obmen_T2M_2(A2C,Pa1,Pc,mu,Cp1,Cp(v_c),r1,RC.ACc,KFA,KFC);
     [A2GL,A2GW]=Obmen_T2M_2(A2G,Pa1,Pg,mu,Cp1,Cp(v_g),r2,RC.AGc,KFA,KFG);
     [A2DL,A2DW]=Obmen_T2M_2(A2D,Pa1,Pd1,mu,Cp1,Cp1d,r3,c3,KFA,KFD);
     
     [D2CL,D2CW]=Obmen_T2M_2(D2C,Pd1,Pc,mu,Cp1d,Cp(v_c),r1d,RC.DCc,KFD,KFC);
     [D2GL,D2GW]=Obmen_T2M_2(D2G,Pd1,Pg,mu,Cp1d,Cp(v_g),r2d,RC.DGc,KFD,KFG);

     [W1,W6,W7]=Well_MKT_2(wom(:,2),wom(:,1),Uf(wom(:,3)),Cp(v_a(v1==1)),mu,CpW(wom(:,3)),kfw(va),kfo(va));
     [W1C,W6C,W7C]=Well_MKT_2(WonC(:,2),WonC(:,1),Uf(WonC(:,3)),Cp(v_c),mu,CpW(WonC(:,3)),kfw(vc),kfo(vc));
     [W1G,W6G,W7G]=Well_MKT_2(WonG(:,2),WonG(:,1),Uf(WonG(:,3)),Cp(v_g),mu,CpW(WonG(:,3)),kfw(vg),kfo(vg));
     [W1D,W6D,W7D]=Well_MKT_2(wod(:,2),wod(:,1),Uf(wod(:,3)),Cp(v_d(v2==1)),mu,CpW(wod(:,3)),kfw(vd),kfo(vd));

     CL1=sparse(RC.Cr2,RC.Cc2,CL,nc,nc);  CL1=CL1+CL1';
     GL1=sparse(RC.Gr2,RC.Gc2,GL,ng,ng);  GL1=GL1+GL1';
     
     CW1=sparse(RC.Cr2,RC.Cc2,CW,nc,nc);  CW1=CW1+CW1';
     GW1=sparse(RC.Gr2,RC.Gc2,GW,ng,ng);  GW1=GW1+GW1';

     A2CL=sparse(r1,RC.ACc,A2CL,na,nc);
     A2GL=sparse(r2,RC.AGc,A2GL,na,ng);
     A2DL=sparse(r3,c3,A2DL,na,nd);
     D2CL=sparse(r1d,RC.DCc,D2CL,nd,nc);
     D2GL=sparse(r2d,RC.DGc,D2GL,nd,ng);
     
     A2CW=sparse(r1,RC.ACc,A2CW,na,nc);
     A2GW=sparse(r2,RC.AGc,A2GW,na,ng);
     A2DW=sparse(r3,c3,A2DW,na,nd);
     D2CW=sparse(r1d,RC.DCc,D2CW,nd,nc);
     D2GW=sparse(r2d,RC.DGc,D2GW,nd,ng);
     
     C1=CL1-sparse(1:nc,1:nc,sum(CL1,1)+sum(A2CL,1)+sum(D2CL,1),nc,nc)-sparse(WonC(:,1),WonC(:,1),W1C,nc,nc);
     G1=GL1-sparse(1:ng,1:ng,sum(GL1,1)+sum(A2GL,1)+sum(D2GL,1),ng,ng)-sparse(WonG(:,1),WonG(:,1),W1G,ng,ng);
     
     C2=CW1-sparse(1:nc,1:nc,sum(CW1,1)+sum(A2CW,1)+sum(D2CW,1),nc,nc)-sparse(WonC(:,1),WonC(:,1),W6C,nc,nc);
     G2=GW1-sparse(1:ng,1:ng,sum(GW1,1)+sum(A2GW,1)+sum(D2GW,1),ng,ng)-sparse(WonG(:,1),WonG(:,1),W6G,ng,ng);

    [bAl,bAw,bl,bw]=Potok_GY(T_gy,Pgy,Pa,rc_gy,KfwM,KfoM,v1,mu,Na);
    [bDl,bDw,bld,bwd]=Potok_GY(T_gy_d,Pgy2,Pd,rc_gy_d,KfwD,KfoD,v2,mu,Nd);
    
     TL2=TL2(:,v1==1);
     TL2=TL2(v1==1,:);

     DL2=DL2(:,v2==1);
     DL2=DL2(v2==1,:);
     
     A1=TL2-sparse(1:na,1:na,sum(TL2,2)+sum(A2CL,2)+sum(A2GL,2)+sum(A2DL,2)+bAl',na,na)-sparse(wom(:,1),wom(:,1),W1,na,na);
     D1=DL2-sparse(1:nd,1:nd,sum(DL2,2)+sum(D2CL,2)+sum(D2GL,2)+sum(A2DL,1)'+bDl',nd,nd)-sparse(wod(:,1),wod(:,1),W1D,nd,nd);
  
     AMC1=[A1,   A2CL,  A2GL,  A2DL;
          A2CL',  C1,   C2GL,  D2CL';
          A2GL', C2GL',  G1,   D2GL';
          A2DL', D2CL,  D2GL,   D1];

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

    WM1=WM1(wom(Qf~=0,3),:);
    WM2=WM2(:,wom(Qf~=0,3));
    WM3=WM3(wom(Qf~=0,3),wom(Qf~=0,3));

    BC1=[ba1-bl';bc1;bg1;bd1-bld']';

     Pt=[BC1,Qz(wom(Qf~=0,3))']/[AMC1,WM2;WM1,WM3];
     Pj(:,1)=Pt(1:na+nc+ng+nd);

     Pw(wom(Qf~=0,3))=Pt(na+nc+ng+nd+1:end);

     TW2=TW2(:,v1==1);
     TW2=TW2(v1==1,:);
          
     DW2=DW2(:,v2==1);
     DW2=DW2(v2==1,:);
     
     A2=TW2-sparse(1:na,1:na,sum(TW2,2)+sum(A2CW,2)+sum(A2GW,2)+sum(A2DW,2)+bAw',na,na)-sparse(wom(:,1),wom(:,1),W6,na,na);
     D2=DW2-sparse(1:nd,1:nd,sum(DW2,2)+sum(D2CW,2)+sum(D2GW,2)+sum(A2DW,1)'+bDw',nd,nd)-sparse(wod(:,1),wod(:,1),W6D,nd,nd);
     
     AMC2=[A2,   A2CW,  A2GW,  A2DW;
          A2CW',  C2,   C2GW,  D2CW';
          A2GW', C2GW',  G2,   D2GW';
          A2DW', D2CW,  D2GW,   D2];
     
     ba2=sparse(wom(:,1),ones(1,size(wom,1)),W6.*Pw(wom(:,3)),na,1);
     bc2=sparse(WonC(:,1),ones(1,size(WonC,1)),W6C.*Pw(WonC(:,3)),nc,1);
     bg2=sparse(WonG(:,1),ones(1,size(WonG,1)),W6G.*Pw(WonG(:,3)),ng,1);
     bd2=sparse(wod(:,1),ones(1,size(wod,1)),W6D.*Pw(wod(:,3)),nd,1);

     ba2=ba2+bw';
     bd2=bd2+bwd';
     BC=[ba2;bc2;bg2;bd2];
     
     [ndt,j_ndt,fl]=vibor_t(ndt,Fl,Pc,Pg,Pw,PR,RC,dt,SW0([v_c,v_g]),SW([v_c,v_g]),CL,GL,dVc,dVg,W1C,W1G,...
         WonC,WonG,i,j_ndt);
    
     Sw1=Sw;
     Sw=Sw1+(AMC2*Pj+BC)./[dVa;dVc;dVg;dVd]*dt/ndt;
  %   Sw(vc)-Sw1(vc)
     Sw=Sw.*(Sw>=0).*(Sw<=1)+(Sw>1);
%fghgfh
     
     SW0([v_c,v_g])=Sw1([vc,vg]);
     SW([v_c,v_g])=Sw([vc,vg]);
     
      %hj(:,i)=SCw(vc);

[Blc,Bwc]=crack_bond(Blc,Bwc,Pj(vc),Pj(va),A2CL,A2CW,dt,ndt,c1);
[Blg,Bwg]=crack_bond(Blg,Bwg,Pj(vg),Pj(va),A2GL,A2GW,dt,ndt,c2);

[Blcd,Bwcd]=crack_bond(Blcd,Bwcd,Pj(vc),Pj(vd),D2CL,D2CW,dt,ndt,c1d);
[Blgd,Bwgd]=crack_bond(Blgd,Bwgd,Pj(vg),Pj(vd),D2GL,D2GW,dt,ndt,c2d);

% full([Pj(won)<Pw(wn),Uf(wn)])
% dfgh
Pc=Pj(vc);
Pg=Pj(vg);

Pa(v1==1)=Pj(va);
Pd(v2==1)=Pj(vd);
Pi=[Pa;Pc;Pg;Pd];

Qm(:,:)=Qm+QBild(W1,W6,W7,Pj(va),Uf(wom(:,3)),wom(:,1),dt/ndt,Pw(wom(:,3)),wom(:,3),nw);
Q1(:,:)=Q1+QBild(W1C,W6C,W7C,Pj(vc),Uf(WonC(:,3)),WonC(:,1),dt/ndt,Pw(WonC(:,3)),WonC(:,3),nw);
Q2(:,:)=Q2+QBild(W1G,W6G,W7G,Pj(vg),Uf(WonG(:,3)),WonG(:,1),dt/ndt,Pw(WonG(:,3)),WonG(:,3),nw);
Qd(:,:)=Qd+QBild(W1D,W6D,W7D,Pj(vd),Uf(wod(:,3)),wod(:,1),dt/ndt,Pw(wod(:,3)),wod(:,3),nw);
%SCwC(Na+1:Na+nc)-1
ndtI(i)=ndt;
% Sw(vc)'
% [Pj(va),Pj(vc(end:-1:1)),Pj(vd)]
end;
%ndtI

%Bc(r1)=Bc;
%hj-1

dSS=sum((Sw([vc,vg])-Sw0([vc,vg])).*[dVc;dVg])+sum(Q1(:,1))+sum(Bwc);
%sum((SCwC(Na+1:end)-Sw0(Na+1:end)).*[dVc;dVg])

CSw=Sw(vc);
GSw=Sw(vg);
end
function [v,won,r1,c1,r2,c2,rc_gy,rc_in_h,T_gy,T_in_h]=ext_cr(CR_rc)
r1=CR_rc.r1;     c1=CR_rc.c1;
r2=CR_rc.r2;     c2=CR_rc.c2;

rc_gy=CR_rc.rc_gy;
rc_in_h=CR_rc.rc_in_h;

T_gy=CR_rc.T_gy;
T_in_h=CR_rc.T_in_h;
won=CR_rc.won;
v=CR_rc.v;
end
function [Bl,Bw]=crack_bond(Bl,Bw,P1,P2,A2CL,A2CW,dt,ndt,c)
 if isempty(c)==0
    Bw=Bw+(A2CW*P1-sum(A2CW,2).*P2)*dt/ndt;
    Bl=Bl+(A2CL*P1-sum(A2CL,2).*P2)*dt/ndt;
 end
end

