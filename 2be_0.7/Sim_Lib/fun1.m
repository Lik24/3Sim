function [Bc,Bg,Blc,Blg,Swc,Swg,i,Q1,Q2,Qm,dSS]=fun1(RC,Pi,SCw,SCp,PR,TC,TG,A2C,A2G,WonC,...
    WonG,Uf,CpW,Pw,dt,dV,CR_rc,Qz,Qf,ndt,Pi0,L,Lc,Lg,Ke)

as=PR.as;
aw=PR.aw;
ts=PR.ts;
tw=PR.tw;
mu=PR.mu;
kms=PR.kms;
Ro=PR.Ro;

na=RC.na;
nc=RC.nc;
ng=RC.ng;

r1=CR_rc(1,1).r1;     c1=CR_rc(1,1).c1;
r2=CR_rc(1,1).r2;     c2=CR_rc(1,1).c2;

rc_gy=CR_rc(1,1).rc_gy;
rc_in_h=CR_rc(1,1).rc_in_h;

T_gy=CR_rc(1,1).T_gy;
%T_in=CR_rc.T_in;
T_in_h=CR_rc(1,1).T_in_h;

r1d=CR_rc(1,2).r1;     c1d=CR_rc(1,2).c1;
r2d=CR_rc(1,2).r2;     c2d=CR_rc(1,2).c2;
rc_gy_d=CR_rc(1,2).rc_gy;
rc_in_hd=CR_rc(1,2).rc_in_h;
TD_in_h=CR_rc(1,2).T_in_h;

won=CR_rc(1,1).won;    wf=CR_rc.wf;  wn1=CR_rc.wn1;   wn=CR_rc.wn;

va=1:na;
vc=na+1:na+nc;
vg=na+nc+1:na+nc+ng;
vd=na+nc+ng+1:na+nc+ng+nd;
 %Pi1=Pi;
 %Pi=Pi0;
Pa=Pi(va);    Pc=Pi(vc);    Pg=Pi(vg);      Pd=Pi(vd);
dVa=dV(va);   dVc=dV(vc);   dVg=dV(vg);     dVd=dV(vd);
 


Qz=Qz/dt;
nw=size(Qz,1);
Qf=Qf(wn1);

%%
vP1=Pa(RC.ACr)>=Pc(RC.ACc);
vP2=Pa(RC.AGr)>=Pg(RC.AGc);
    fl1=sum([SCw(vc);SCw(vg)])/size([SCw(vc);SCw(vg)],1)==1;
    fl2=sum([vP1;vP2])==0;
    Fl=fl1*fl2;
    
%%
C2G=sparse(nc,ng);
C2GW=C2G;

Bc=zeros(size(A2C,2),1);
Bg=zeros(size(A2G,2),1);
Blc=zeros(size(A2C,2),1);
Blg=zeros(size(A2G,2),1);

Q1=zeros(nw,5);
Q2=zeros(nw,5);
Qm=zeros(nw,5);

 v1=zeros(na,1);
 v1([RC.ACr;RC.AGr])=1;
 Na=sum(v1);
 dVa(v1~=1)=[];

 v2=zeros(nd,1);
 v2([RC.DCr;RC.DGr])=1;
 Nd=sum(v2);
 dVd(v2~=1)=[];
 
va1=1:Na;
vc1=Na+1:na+nc;
vg1=Na+nc+1:na+nc+ng;
vd1=Na+nc+ng+1:na+nc+ng+Nd;
 
Pj(:,1)=[Pi(v1==1);Pc;Pg;Pd(v2==1);];
Pgy=Pa(rc_gy(:,1));
Pgy2=Pd(rc_gy_d(:,1));
 
Swa=SCw(va);
Swd=SCw(vd);
SCwC=[Swa(v1==1);SCw(vc);SCw(vg);Swd(v2==1)];
Sw0=SCwC;
SCw0=SCw;

vac=sort(r1);
vca=Na+sort(c1);
vag=sort(r2);
vga=Na+nc+sort(c2);
% vca
% Na
kfw=zeros(size(SCwC));
kfo=zeros(size(SCwC));

Kfw=Sat_cal(Swa,1,1,as,aw); %water
Kfo=Sat_cal(Swa,2,1,as,aw); %oil

fl=1;
i=0;
j_ndt=0;

while fl==1% & i<3000
i=i+1;

     kfw(1:Na)=Sat_cal(SCwC(1:Na),1,1,as,aw); %water
     kfo(1:Na)=Sat_cal(SCwC(1:Na),2,1,as,aw); %oil
     
     kfw(Na+1:end)=Sat_cal(SCwC(Na+1:end),1,1,ts,tw); %water
     kfo(Na+1:end)=Sat_cal(SCwC(Na+1:end),2,1,ts,tw); %oil
     
     Kfw(v1==1)=kfw(1:Na);
     Kfo(v1==1)=kfo(1:Na);
     
     [vPa1,vPc1,vPc2,vPg1,vPg2,vPd1,dPa,dPc,dPg,dPd]=pre_potok_2(Pi,Pj,RC,rc_in_h,rc_in_hd,Na,Nd);
     
     %[TL2,TW2]=Potok_MKT_2(T_in,vPa1,SCp(va),mu,rc_in,na,Kfw,Kfo);
     [TL2,TW2]=Potok_MKT_2(T_in_h,vPa1,SCp(va),mu,rc_in_h,na,Kfw,Kfo,kms(1),dPa,L,Ro,Ke);
     [CL,CW]=Potok_Tube_2(TC,Pc,vPc1,vPc2,kfw(vca),kfo(vca),SCp(vc),PR,RC.Cr2,RC.Cc2,kms(2),dPc,Lc,nc);
     [GL,GW]=Potok_Tube_2(TG,Pg,vPg1,vPg2,kfw(vga),kfo(vga),SCp(vg),PR,RC.Gr2,RC.Gc2,kms(3),dPg,Lg,ng);
     [DL2,DW2]=Potok_MKT_2(TD_in_h,vPd1,SCp(vd),mu,rc_in_h,nd,kfw(vda),kfo(vda),kms(4),dPd,L,Ro,Ke);

     Pa11=Pi(va(v1==1));
     Cp11=SCp(va(v1==1));

     KFA=[kfw(1:Na),kfo(1:Na)];
     KFC=[kfw(vca),kfo(vca)];
     KFG=[kfw(vga),kfo(vga)];
     KFD=[kfw(vda),kfo(vda)];
     
     [A2CL,A2CW]=Obmen_T2M_2(A2C,Pa11,Pc,mu,Cp11,SCp(vc),r1,RC.ACc,KFA,KFC);
     [A2GL,A2GW]=Obmen_T2M_2(A2G,Pa11,Pg,mu,Cp11,SCp(vg),r2,RC.AGc,KFA,KFG);
     [A2DL,A2DW]=Obmen_T2M_2(A2D,Pa11,Pd,mu,Cp11,SCp(vd),r2d,RC.ADc,KFA,KFD);
     
     [D2CL,D2CW]=Obmen_T2M_2(D2C,Pd11,Pc,mu,Cp1d,SCp(vc),r1d,RC.DCc,KFD,KFC);
     [D2GL,D2GW]=Obmen_T2M_2(D2G,Pd11,Pg,mu,Cp1d,SCp(vg),r2d,RC.DGc,KFD,KFG);

     [W1,W6,W7]=Well_MKT_2(wf,won,Uf(wn),SCp(va(v1==1)),mu,CpW(wn),kfw(1:Na),kfo(1:Na));
     [W1C,W6C,W7C]=Well_MKT_2(WonC(:,2),WonC(:,1),Uf(WonC(:,3)),SCp(vc),mu,CpW(WonC(:,3)),kfw(vca),kfo(vca));
     [W1G,W6G,W7G]=Well_MKT_2(WonG(:,2),WonG(:,1),Uf(WonG(:,3)),SCp(vg),mu,CpW(WonG(:,3)),kfw(vga),kfo(vga));
     [W1D,W6D,W7D]=Well_MKT_2(WonD(:,2),WonD(:,1),Uf(WonD(:,3)),SCp(vd),mu,CpW(WonD(:,3)),kfw(vda),kfo(vda));

     CL1=sparse(RC.Cr2,RC.Cc2,CL,nc,nc);  CL1=CL1+CL1';
     GL1=sparse(RC.Gr2,RC.Gc2,GL,ng,ng);  GL1=GL1+GL1';
     
     CW1=sparse(RC.Cr2,RC.Cc2,CW,nc,nc);  CW1=CW1+CW1';
     GW1=sparse(RC.Gr2,RC.Gc2,GW,ng,ng);  GW1=GW1+GW1';

     A2CL=sparse(r1,RC.ACc,A2CL,Na,nc);
     A2GL=sparse(r2,RC.AGc,A2GL,Na,ng);
     
     A2CW=sparse(r1,RC.ACc,A2CW,Na,nc);
     A2GW=sparse(r2,RC.AGc,A2GW,Na,ng);

     C1=CL1-sparse(1:nc,1:nc,sum(CL1,1)+sum(A2CL,1),nc,nc)-sparse(WonC(:,1),WonC(:,1),W1C,nc,nc);
     G1=GL1-sparse(1:ng,1:ng,sum(GL1,1)+sum(A2GL,1),ng,ng)-sparse(WonG(:,1),WonG(:,1),W1G,ng,ng);
     
     C2=CW1-sparse(1:nc,1:nc,sum(CW1,1)+sum(A2CW,1),nc,nc)-sparse(WonC(:,1),WonC(:,1),W6C,nc,nc);
     G2=GW1-sparse(1:ng,1:ng,sum(GW1,1)+sum(A2GW,1),ng,ng)-sparse(WonG(:,1),WonG(:,1),W6G,ng,ng);

    [bAl,bAw,bl,bw]=Potok_GY(T_gy,Pgy,Pa,rc_gy,Kfw,Kfo,v1,mu,na);
    [bDl,bDw,bld,bwd]=Potok_GY(T_gy_d,Pgy,Pd,rc_gy_d,Kfw,Kfo,v2,mu,nd);
    
     TL2=TL2(:,v1==1);
     TL2=TL2(v1==1,:);

     DL2=DL2(:,v2==1);
     DL2=DL2(v2==1,:);
     
     A1=TL2-sparse(1:Na,1:Na,sum(TL2,2)+sum(A2CL,2)+sum(A2GL,2)+sum(A2DL,2)+bAl',Na,Na)-sparse(won,won,W1,Na,Na);
     D1=DL2-sparse(1:Nd,1:Nd,sum(DL2,2)+sum(D2CL,2)+sum(D2GL,2)+sum(A2DL,2)+bDl',Nd,Nd)-sparse(won,won,W1,Nd,Nd);
  
     AMC1=[A1,   A2CL,  A2GL,   A2DL;
          A2CL',  C1,   C2GL,   C2DL;
          A2GL', C2GL',  G1,    G2DL;
          A2DL', C2DL', G2DL',   D1];

    PwNl=Pw(wn);
    ba1=sparse(won,ones(1,size(won,1)),-W1.*PwNl,Na,1);
    bc1=sparse(WonC(:,1),ones(1,size(WonC,1)),-W1C.*Pw(WonC(:,3)),nc,1);
    bg1=sparse(WonG(:,1),ones(1,size(WonG,1)),-W1G.*Pw(WonG(:,3)),ng,1);
    bd1=sparse(WonD(:,1),ones(1,size(WonD,1)),-W1D.*Pw(WonD(:,3)),nd,1);

    W2M=sparse(wn,won,W1,nw,Na);
    W2C=sparse(WonC(:,3),WonC(:,1),W1C,nw,nc);
    W2G=sparse(WonG(:,3),WonG(:,1),W1G,nw,ng);
    W2D=sparse(WonD(:,3),WonD(:,1),W1D,nw,nd);
    
    ba1=ba1.*(sum(W2M(Qz==0,:),1)~=0)';
    bc1=bc1.*(sum(W2C(Qz==0,:),1)~=0)';
    bg1=bg1.*(sum(W2G(Qz==0,:),1)~=0)';
    bd1=bd1.*(sum(W2D(Qz==0,:),1)~=0)';
    
    WM1=[W2M,W2C,W2G,W2D];
    WM2=WM1';
    W3vec=sparse(wn,1,W1,nw,1)+sparse(WonC(:,3),1,W1C,nw,1)+sparse(WonG(:,3),1,W1G,nw,1)+sparse(WonD(:,3),1,W1D,nw,1);
    WM3=-sparse(1:nw,1:nw,W3vec,nw,nw); 
    %WM1(Qz~=0,:)
   % full([Qz(wn(Qf~=0)),W3vec(wn(Qf~=0)),Pw(wn(Qf~=0))])

    WM1=WM1(wn1(Qf~=0),:);
    WM2=WM2(:,wn1(Qf~=0));
    WM3=WM3(wn1(Qf~=0),wn1(Qf~=0));

    BC1=[ba1-bl';bc1;bg1;bd1-bld]';

     Pt=[BC1,Qz(wn1(Qf~=0))']/[AMC1,WM2;WM1,WM3];
     Pj(:,1)=Pt(1:Na+nc+ng);

     Pw(wn1(Qf~=0))=Pt(Na+nc+ng+1:end);

     TW2=TW2(:,v1==1);
     TW2=TW2(v1==1,:);
     
     A2=TW2-sparse(1:Na,1:Na,sum(TW2,2)+sum(A2CW,2)+sum(A2GW,2)+sum(A2DW,2)+bAw',Na,Na)-sparse(won,won,W6,Na,Na);
     D2=DW2-sparse(1:Na,1:Na,sum(DW2,2)+sum(D2CW,2)+sum(D2GW,2)+sum(A2DW,2)+bDw',Nd,Nd)-sparse(won,won,W6D,Nd,Nd);
     
     AMC2=[A2,   A2CW,  A2GW,   A2DW;
          A2CW',  C2,   C2GW,   C2DW;
          A2GW', C2GW',  G2,    G2DW;
          A2DW', C2DW', G2DW',   D2];
     
     ba2=sparse(won(:,1),ones(1,size(won,1)),W6.*Pw(wn),Na,1);
     bc2=sparse(WonC(:,1),ones(1,size(WonC,1)),W6C.*Pw(WonC(:,3)),nc,1);
     bg2=sparse(WonG(:,1),ones(1,size(WonG,1)),W6G.*Pw(WonG(:,3)),ng,1);
     bd2=sparse(WonD(:,1),ones(1,size(WonD,1)),W6D.*Pw(WonD(:,3)),nd,1);

     ba2=ba2+bw';
     bd2=bd2+bwd';
     BC=[ba2;bc2;bg2;bd2];
%      fCL=full(CL);
%      fGL=full(GL);
     
%       [ndt,j_ndt,fl]=vibor_t_mex(ndt,Fl,full(Pc),full(Pg),Pw,PR,RC,dt,SCw0([vc,vg]),SCw([vc,vg]),fCL,fGL,dVc,dVg,W1C,W1G,...
%           WonC,WonG,i,j_ndt);
     
     [ndt,j_ndt,fl]=vibor_t(ndt,Fl,Pc,Pg,Pw,PR,RC,dt,SCw0([vc,vg]),SCw([vc,vg]),CL,GL,dVc,dVg,W1C,W1G,...
         WonC,WonG,i,j_ndt);
    
     SCwC1=SCwC;
     SCwC=SCwC1+(AMC2*Pj+BC)./[dVa;dVc;dVg;dVd]*dt/ndt;
     SCwC=SCwC.*(SCwC>=0).*(SCwC<=1)+(SCwC>1);
%fghgfh
     
     SCw0(na+1:end)=SCwC1(Na+1:end);
     SCw(na+1:end)=SCwC(Na+1:end);
     
      %hj(:,i)=SCw(vc);

Bc=Bc+(A2CW*Pj(vc1)-sum(A2CW,2).*Pj(va1)+D2CW*Pj(vd1)-sum(D2CW,2).*Pj(vd))*dt/ndt;
Blc=Blc+(A2CL*Pj(Na+1:Na+nc)-sum(A2CL,2).*Pj(1:Na)+D2CL*Pj(Na+nc+ng+1:Na+nc+ng+Nd)-sum(D2CL,2).*Pj(vd))*dt/ndt;

if isempty(c2)==0
    Bg=Bg+(A2GW*Pj(Na+nc+1:end)-sum(A2GW,2).*Pj(1:Na)+D2GW*Pj(Na+nc+1:end)-sum(D2GW,2).*Pj(vd))*dt/ndt;
    Blg=Blg+(A2GL*Pj(Na+nc+1:end)-sum(A2GL,2).*Pj(1:Na))*dt/ndt;
end;

% full([Pj(won)<Pw(wn),Uf(wn)])
% dfgh
Pc=Pj(Na+1:Na+nc);
Pg=Pj(Na+nc+1:Na+nc+ng);

Pa(v1==1)=Pj(1:Na);
Pi=[Pa;Pc;Pg];

Qm(:,:)=Qm+QBild(W1,W6,W7,Pj(1:Na),Uf(wn),won,dt/ndt,Pw(wn),wn,nw);
Q1(:,:)=Q1+QBild(W1C,W6C,W7C,Pj(Na+1:Na+nc,1),Uf(WonC(:,3)),WonC(:,1),dt/ndt,Pw(WonC(:,3)),WonC(:,1),nw);
Q2(:,:)=Q2+QBild(W1G,W6G,W7G,Pj(vga,1),Uf(WonG(:,3)),WonG(:,1),dt/ndt,Pw(WonG(:,3)),WonG(:,1),nw);
Qd(:,:)=Qd+QBild(W1D,W6D,W7D,Pj(Na+nc+ng+1:Nd),Uf(wn),WonD(:,1),dt/ndt,Pw(wn),wn,nw);
%SCwC(Na+1:Na+nc)-1
ndtI(i)=ndt;
end;
%ndtI

%Bc(r1)=Bc;
%hj-1

dSS=sum((SCwC(Na+1:end)-Sw0(Na+1:end)).*[dVc;dVg])+sum(Q1(:,1))+sum(Bc);
%sum((SCwC(Na+1:end)-Sw0(Na+1:end)).*[dVc;dVg])

Swc=SCwC(Na+1:Na+nc);
Swg=SCwC(Na+nc+1:Na+nc+ng);



