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

r1=CR_rc.r1;     c1=CR_rc.c1;
r2=CR_rc.r2;     c2=CR_rc.c2;
%sd=CR_rc.sd;
rc_gy=CR_rc.rc_gy;
%rc_in=CR_rc.rc_in;
rc_in_h=CR_rc.rc_in_h;

T_gy=CR_rc.T_gy;
%T_in=CR_rc.T_in;
T_in_h=CR_rc.T_in_h;
TD_in_h=CR_rc.D_in_h;

won=CR_rc.won;    wf=CR_rc.wf;  wn1=CR_rc.wn1;   wn=CR_rc.wn;

va=1:na;
vc=na+1:na+nc;
vg=na+nc+1:na+nc+ng;
 %Pi1=Pi;
 %Pi=Pi0;
 Pc=Pi(vc);  Pg=Pi(vg);
 dVa=dV(va);  dVc=dV(vc);  dVg=dV(vg);
 
Pa=Pi(va);

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
 %    sort(RC.ACr)
 Na=sum(v1);
 dVa(v1~=1)=[];

Pj(:,1)=[Pi(v1==1);Pc;Pg];
Pgy=Pa(rc_gy(:,1));
 
Swa=SCw(va);
SCwC=[Swa(v1==1);SCw(vc);SCw(vg)];
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
     
     [vPa1,vPc1,vPc2,vPg1,vPg2,dPa,dPc,dPg]=pre_potok_2(Pi,Pj,RC,rc_in_h,Na);
     
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
     
     [A2CL,A2CW]=Obmen_T2M_2(A2C,Pa11,Pc,mu,Cp11,SCp(vc),r1,RC.ACc,KFA,KFC);
     [A2GL,A2GW]=Obmen_T2M_2(A2G,Pa11,Pg,mu,Cp11,SCp(vg),r2,RC.AGc,KFA,KFG);

     [W1,W6,W7]=Well_MKT_2(wf,won,Uf(wn),SCp(va(v1==1)),mu,CpW(wn),kfw(1:Na),kfo(1:Na));
     [W1C,W6C,W7C]=Well_MKT_2(WonC(:,2),WonC(:,1),Uf(WonC(:,3)),SCp(na+1:na+nc),mu,CpW(WonC(:,3)),kfw(vca),kfo(vca));
     [W1G,W6G,W7G]=Well_MKT_2(WonG(:,2),WonG(:,1),Uf(WonG(:,3)),SCp(na+nc+1:end),mu,CpW(WonG(:,3)),kfw(vga),kfo(vga));

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
    
     TL2=TL2(:,v1==1);
     TL2=TL2(v1==1,:);

     A1=TL2-sparse(1:Na,1:Na,sum(TL2,2)+sum(A2CL,2)+sum(A2GL,2)+bAl',Na,Na)-sparse(won,won,W1,Na,Na);
  
     AMC1=[A1,A2CL,A2GL;A2CL',C1,C2GW;A2GL',C2GW',G1];

    PwNl=Pw(wn);
    ba1=sparse(won,ones(1,size(won,1)),-W1.*PwNl,Na,1);
    bc1=sparse(WonC(:,1),ones(1,size(WonC,1)),-W1C.*Pw(WonC(:,3)),nc,1);
    bg1=sparse(WonG(:,1),ones(1,size(WonG,1)),-W1G.*Pw(WonG(:,3)),ng,1);

    W2M=sparse(wn,won,W1,nw,Na);
    W2C=sparse(WonC(:,3),WonC(:,1),W1C,nw,nc);
    W2G=sparse(WonG(:,3),WonG(:,1),W1G,nw,ng);
    
    ba1=ba1.*(sum(W2M(Qz==0,:),1)~=0)';
    bc1=bc1.*(sum(W2C(Qz==0,:),1)~=0)';
    bg1=bg1.*(sum(W2G(Qz==0,:),1)~=0)';
    
    WM1=[W2M,W2C,W2G];
    WM2=WM1';
    W3vec=sparse(wn,1,W1,nw,1)+sparse(WonC(:,3),1,W1C,nw,1)+sparse(WonG(:,3),1,W1G,nw,1);
    WM3=-sparse(1:nw,1:nw,W3vec,nw,nw); 
    %WM1(Qz~=0,:)
   % full([Qz(wn(Qf~=0)),W3vec(wn(Qf~=0)),Pw(wn(Qf~=0))])

    WM1=WM1(wn1(Qf~=0),:);
    WM2=WM2(:,wn1(Qf~=0));
    WM3=WM3(wn1(Qf~=0),wn1(Qf~=0));

    BC1=[ba1-bl';bc1;bg1]';

     Pt=[BC1,Qz(wn1(Qf~=0))']/[AMC1,WM2;WM1,WM3];
     Pj(:,1)=Pt(1:Na+nc+ng);

     Pw(wn1(Qf~=0))=Pt(Na+nc+ng+1:end);

     TW2=TW2(:,v1==1);
     TW2=TW2(v1==1,:);
     A2=TW2-sparse(1:Na,1:Na,sum(TW2,2)+sum(A2CW,2)+sum(A2GW,2)+bAw',Na,Na)-sparse(won,won,W6,Na,Na);
     
     AMC2=[A2,A2CW,A2GW;A2CW',C2,C2GW;A2GW',C2GW',G2];
     
     ba2=sparse(won,ones(1,size(won,1)),W6.*Pw(wn),Na,1);
     bc2=sparse(WonC(:,1),ones(1,size(WonC,1)),W6C.*Pw(WonC(:,3)),nc,1);
     bg2=sparse(WonG(:,1),ones(1,size(WonG,1)),W6G.*Pw(WonG(:,3)),ng,1);

     ba2=ba2+bw';
     BC=[ba2;bc2;bg2];
%      fCL=full(CL);
%      fGL=full(GL);
     
%       [ndt,j_ndt,fl]=vibor_t_mex(ndt,Fl,full(Pc),full(Pg),Pw,PR,RC,dt,SCw0([vc,vg]),SCw([vc,vg]),fCL,fGL,dVc,dVg,W1C,W1G,...
%           WonC,WonG,i,j_ndt);
     
     [ndt,j_ndt,fl]=vibor_t(ndt,Fl,Pc,Pg,Pw,PR,RC,dt,SCw0([vc,vg]),SCw([vc,vg]),CL,GL,dVc,dVg,W1C,W1G,...
         WonC,WonG,i,j_ndt);
    
     SCwC1=SCwC;
     SCwC=SCwC1+(AMC2*Pj+BC)./[dVa;dVc;dVg]*dt/ndt;
     
     SCwC=SCwC.*(SCwC>=0).*(SCwC<=1)+(SCwC>1);
%fghgfh
     
     SCw0(na+1:end)=SCwC1(Na+1:end);
     SCw(na+1:end)=SCwC(Na+1:end);
     
      %hj(:,i)=SCw(vc);

Bc=Bc+(A2CW*Pj(Na+1:Na+nc)-sum(A2CW,2).*Pj(1:Na))*dt/ndt;
Blc=Blc+(A2CL*Pj(Na+1:Na+nc)-sum(A2CL,2).*Pj(1:Na))*dt/ndt;

if isempty(c2)==0
    Bg=Bg+(A2GW*Pj(Na+nc+1:end)-sum(A2GW,2).*Pj(1:Na))*dt/ndt;
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



