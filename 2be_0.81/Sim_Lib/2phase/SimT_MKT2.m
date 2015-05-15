function [Pj,Swj,Soj,Tj,MCpj,p,Q,Pw,PpW,CSwj,NDT,Uf,dt1,dV0,ka,dtz,DSwj]=SimT_MKT2(PR,C,A2C,GData,BB,A2B,DData,dVC,dVB,DATA,WData,GYData,fll,CR_GRUP)
tic

KX=DATA.gKX;
KY=DATA.gKY;
KZ=DATA.gKZ;
Mp=DATA.gMp;
P=DATA.gP;
MSw=DATA.gSw;
MSo=DATA.gSo;    
MCp=DATA.gCp;
XY=DATA.XY;
BND=DATA.BND;
BndXY=DATA.BndXY;
BndZ=DATA.BndZ;
H=DATA.gH;
Z=DATA.gZ;
ka=DATA.ka;

[G,M2FR.A2G,dVG,WELL.WonG,GEOM.Lg,Mp_g]=Rasp(GData,3);   %���������� ��������� ��� ���. ����.
[D,M2FR.A2D,dVD,WELL.WonD,GEOM.Ld,Mp_d]=Rasp(DData,4);   %���������� ��������� ��� ������� �����

WELL.Won=DATA.Won;       %������ ����� �������, 
WELL.WonC=DATA.WonV;

Pw=WData.Pw;
WELL.Uf=WData.Uf;
WELL.CpW=WData.CpW;
TeW=WData.TeW;
Qz=WData.Qz;
PwQC_bnd=WData.PwQC_bnd;

Nl=PR.Nl;  as=PR.as; aw=PR.aw; ts=PR.ts; tw=PR.tw; mu=PR.mu; mup=PR.mup;
Ta=PR.Ta;  dtt=PR.dt; ndt=PR.ndt; zc=PR.zc; Bo=PR.Bo;  kms=PR.kms;  Ro=PR.Ro;
rs=PR.rs; nw=size(Qz,1);

NCELL.nw = nw;

Qm=zeros(nw,5,size(WELL.Uf,2));
Qc=zeros(nw,5,size(WELL.Uf,2));
Qg=zeros(nw,5,size(WELL.Uf,2));
Qd=zeros(nw,5,size(WELL.Uf,2));

sQm1=zeros(size(WELL.Uf,1),5);
sQ1=zeros(size(WELL.Uf,1),5);
sQ2=zeros(size(WELL.Uf,1),5);
sQd1=zeros(size(WELL.Uf,1),5);

[A]=MR_Prop_Bond(XY,Nl,BND);
[GEOM.L,B,S,H1,HV]=Geome3_1(A,XY,Z,H);

[r,c]=find(A(1:size(XY,1),1:size(XY,1))==1);

Wf=KWell(KX,H,S,GEOM.L,B,WELL.Won,r,c,WData.Doly,WData.SDoly,WData.r0,XY,nw,Nl);
WELL.Won(:,2)=KWell_Horiz(Wf,KX,KY,KZ,H,S,GEOM.L,B,WELL.Won,r,c,WData.Doly,WData.SDoly,WData.r0,XY,nw,Nl);

XY=repmat(XY,Nl,1);
ka1=ka(WELL.Won(:,1));

[A,GEOM.L,B,S,H1,HV,XY,KX,KY,KZ,Mp,MSw,MSo,H,Z,P,MCp,WELL.Won]=YdalActcell(A,GEOM.L,B,S,H1,HV,XY,KX,KY,KZ,Mp,MSw,MSo,H,Z,P,MCp,WELL.Won,ka);
[A,GEOM.L,S,B,H1,HV,GEOM.K,XY,Mp,MSw,MSo,H,Z,P,MCp,T,NTG,p,rz,cz,BXY,BZ,dH,NL,NamXY,GYData]=PereYpor(A,GEOM.L,S,B,H1,HV,KX,KY,KZ,Mp,MSw,MSo,...
    XY,H,Z,P,MCp,DATA,GYData,ka);

[r,c]=find(GEOM.L);
nw1=size(p,2);
nwc=size(A2C,1);
raz=nw1-nwc;

M2FR.A2C = A2C;
M2FR.A2B = A2B;

M2FR.A2C(raz+1:end+raz,:)=M2FR.A2C;  %A2C - ����� ��� � ������������� ���������
M2FR.A2C(1:raz,:)=0;
M2FR.A2C=M2FR.A2C(p,:);
M2FR.A2G=M2FR.A2G(p,:);  % - � ���������������
M2FR.A2D=M2FR.A2D(p,:);  % - ������ ����� (������ �������)
M2FR.A2B=M2FR.A2B(p,:);

for i=1:size(WELL.Won,1)
    WELL.Won(i,1)=find(WELL.Won(i,1)==p);
end;

na=size(XY,1);  nc=size(C,1);  ng=size(G,1);  nd=size(D,1); nb=size(BB,1); nw=size(Qz,1);
% na - ����� ����� � ������� �������, nc - ����� ����� � ��������, ng - ����� ����� � �����. ��������, 
% nd - ����� ����� � ������� �����, nb - ����� ����� � ��������

va=1:na;                                VEC.va=va; % - ������� �������
vc=na+1:na+nc;                          VEC.vc=vc; % - �������
vg=na+nc+1:na+nc+ng;                    VEC.vg=vg; % - �������������� �������
vd=na+nc+ng+1:na+nc+ng+nd;              VEC.vd=vd; % - ������� �����
vb=na+nc+ng+nd+1:na+nc+ng+nd+nb;        VEC.vb=vb; % - �������
VEC.v = [va,vc,vg,vd];

M2FR.C2G=sparse(nc,ng);     M2FR.C2GL=M2FR.C2G;
M2FR.C2B=sparse(nc,nb);     M2FR.C2BL=M2FR.C2B;
M2FR.G2B=sparse(ng,nb);     M2FR.G2BL=M2FR.G2B;
M2FR.D2B=sparse(nd,nb);
b1wb=sparse(nb,1);

M2FR.C2D=DData.C2D;     M2FR.D2CL=M2FR.C2D;
M2FR.G2D=DData.G2D;     M2FR.D2GL=M2FR.G2D;

Swr(va,1)=PR.aw(4); Swr([vc,vg,vd],1)=PR.tw(4);
Sor(va,1)=PR.aw(5); Sor([vc,vg,vd],1)=PR.tw(5);
PR.Swr=Swr; PR.Sor=Sor;

[GEOM.Ke,Ke_gy,dV]=KH2Mat(GEOM.K,HV,Mp,S,r,c,rz,cz,GYData.GY_Kz,GYData.GY_Kxy,BndXY(p),DData);

Mp_c=ones(nc,1);
Mp_b=100*ones(nb,1);

CCp(:,1)=zeros(nc,1);
GCp(:,1)=zeros(ng,1);
DCp(:,1)=zeros(nd,1);
BCp(:,1)=zeros(nb,1);

dVCG=[dV;dVC;dVG;dVD.*Mp_d;dVB];
dV0=[dV;dVC;dVG;dVD.*Mp_d];

[TRM,RC,dV1,dV2]=Pre_fast(A,C,G,D,M2FR.A2C,M2FR.A2G,M2FR.A2D,M2FR.C2G,M2FR.C2D,M2FR.G2D,...
    GEOM.Ke,GEOM.L,B,S,H1,GEOM.K(:,1),DData.Kd,Ke_gy,BndXY(p),BndZ(p),nb,dVCG);  %������������ �������������� �������������
[dZ]=dZ_bild(Z,RC,PR.g,PR.Ro);

vad=RC.ADr;
CSw(:,1)=MSw(RC.ACr,1);
GSw(:,1)=MSw(RC.AGr,1);
DSw(:,1)=MSw(RC.ADr,1);
BSw(:,1)=ones(nb,1);

CSo(:,1)=MSo(RC.ACr,1);
GSo(:,1)=MSo(RC.AGr,1);
DSo(:,1)=MSo(RC.ADr,1);     %///////////////
BSo(:,1)=ones(nb,1);

Pi(va,1)=P;
Pi(vc,1)=P(RC.ACr);
Pi(vg,1)=P(RC.AGr);
Pi(vd,1)=P(RC.ADc);
Pi(vb,1)=GYData.P0(vb);

Ti(va,1)=T;
Ti(vc,1)=T(RC.ACr);
Ti(vg,1)=T(RC.AGr);
Ti(vd,1)=T(RC.ADr);
%Ti(vb,1)=GYData.T0;

Phi(VEC.va,:) = [Pi(VEC.va) - Z(VEC.va)*Ro(1)*PR.g, Pi(VEC.va) - Z(VEC.va)*Ro(2)*PR.g];
Phi(VEC.vc,:) = [Pi(VEC.vc) - Z(RC.ACr)*Ro(1)*PR.g, Pi(VEC.vc) - Z(RC.ACr)*Ro(2)*PR.g];
Phi(VEC.vg,:) = [Pi(VEC.vg) - Z(RC.AGr)*Ro(1)*PR.g, Pi(VEC.vg) - Z(RC.AGr)*Ro(2)*PR.g];
Phi(VEC.vd,:) = [Pi(VEC.vd) - Z(RC.ADr)*Ro(1)*PR.g, Pi(VEC.vd) - Z(RC.ADr)*Ro(2)*PR.g];
Phi(VEC.vb,:) = [Pi(VEC.vb), zeros(nb,1), zeros(nb,1)];
    
j=0;
Pwt=zeros(size(Pw));
PpW=zeros(size(WELL.Uf));

% [CR_rc]=Pre_Crack1(RC,na,TM,A2C,A2G,Wf,Won,WonM,nw);
[CR_rc]=Pre_Crack(RC,na,nd,TRM.TTM,TRM.TD,M2FR.A2C,M2FR.A2G,M2FR.C2D,M2FR.G2D,M2FR.A2D,WELL.Won,WELL.WonD,dZ);

%[CR_ind]=Pre_Crack_p(RC,na,TM,Wf,Won,WonM,nw,CR_GRUP,C,G,A2C,A2G,K(:,1),WonC,WonG);
[CR]=SS_ind(RC,na);

st=0;
t=0;
t_flag=1;
dt=dtt;

CMP.Bo0 = Bo;
CMP.Mp0 = [Mp.*ones(na,1);Mp_c;Mp_g;Mp_d;Mp_b];

GEOM.dV = dVCG./CMP.Mp0; 

Sw=[MSw(:,1);CSw(:,1);GSw(:,1);DSw(:,1);BSw(:,1)];
So=[MSo(:,1);CSo(:,1);GSo(:,1);DSo(:,1);BSo(:,1)];  %//////////
Cp=[MCp(:,1);CCp(:,1);GCp(:,1);DCp(:,1);BCp(:,1)];

V0=sum(dVCG.*(1-Sw));
CMP.P0=Pi;

dtt1=dtt;
Wf0=WELL.Won(:,2);
zapt=1;
 
dt = 1; 
dt1(1) = dt;
Pwt(:,1) = Pw(:,1);
[CMP]=SGim20(Sw,So,zc,Pi,CMP,GEOM.dV,dt);
        
 BXYZ.mxy = find(TRM.Txyz_GY_A(:,1));
 BXYZ.dxy = find(TRM.Txyz_GY_D(:,1));
 BXYZ.mz = find(TRM.Txyz_GY_A(:,2));
 BXYZ.dz = find(TRM.Txyz_GY_A(:,2)); 
 ndt = 100; 
 TL = ones(na,na);
 DL = ones(nd,nd);
while t_flag==1
    t=t+1;
%    if t==1 && dtt==0
%        dt=1;
%        dtt=1;
%    else
%        dtt=dtt1;
%    end;
    
    ft=floor(st) + 1;
    % if ~isempty(find(CpW(:,t)~=0)) fp=1; else fp=0; end;
    fp=0;
    Qf=Qz(:,ft);
    
    Sw0 = Sw;
    So0 = So;
    CMP.Bw(:,1) = CMP.Bw(:,2);
    CMP.Bo(:,1) = CMP.Bo(:,2);
    CMP.Mp(:,1) = CMP.Mp(:,2);
    
   [KWOG.w,KWOG.o] = Sat_cal(Sw,1,as,aw); %water
    
   [W1,~,~,~,QQ.QQm,QQ.QQmwo]=Well_MKT20(Pi(VEC.va),Pwt(WELL.Won(:,3),t),WELL.Won,WELL.Uf(WELL.Won(:,3)),Cp(VEC.va),PR,WELL.CpW(WELL.Won(:,3)),CMP,KWOG,VEC.va,na);% W1 - ������������ �� ���� ��������, W6 - ������ ��� ����, W7 - ������� 
   [W1C,~,~,~,QQ.QQc,QQ.QQcwo]=Well_MKT20(Pi(VEC.vc),Pwt(WELL.WonC(:,3),t),WELL.WonC,WELL.Uf(WELL.WonC(:,3)),Cp(VEC.vc),PR,WELL.CpW(WELL.WonC(:,3)),CMP,KWOG,VEC.vc,nc);
   [W1G,~,~,~,QQ.QQg,QQ.QQgwo]=Well_MKT20(Pi(VEC.vg),Pwt(WELL.WonG(:,3),t),WELL.WonG,WELL.Uf(WELL.WonG(:,3)),Cp(VEC.vg),PR,WELL.CpW(WELL.WonG(:,3)),CMP,KWOG,VEC.vg,ng);
   [W1D,~,~,~,QQ.QQd,QQ.QQdwo]=Well_MKT20(Pi(VEC.vd),Pwt(WELL.WonD(:,3),t),WELL.WonD,WELL.Uf(WELL.WonD(:,3)),Cp(VEC.vd),PR,WELL.CpW(WELL.WonD(:,3)),CMP,KWOG,VEC.vd,nd);
   
   [WBND,QQBND,QQoBND,KWOG_GY]=GY_bild20(GYData,Pi(:,1),Sw(:,1),Cp(:,1),RC,TRM.Txyz_GY_A,TRM.Txyz_GY_D,PR,CMP,BXYZ,VEC);
    dt=vibor_t2(dtt,Pi,RC,dVCG,TL,W1,WELL,Pwt(:,t),na,PR,st,Ta,Sw,Sw0,dt,Nl,va,vd,DL,W1D,nd,dV1,dV2);   
 
    if isempty(RC.Cr)==0 || isempty(RC.Gr)==0
      [Pi,Phi,Sw,Sw0,Pw(:,ft),CMP,ndt,QQ,sQc,sQg]=fun2(RC,Pi,Phi,Sw,Cp,PR,TRM,M2FR,WELL,Pwt(:,t),dt,CR_rc,Qf,ndt,GEOM,DATA.Lc,CMP,ft,QQ,KWOG);
        
   %    [Pi(va),Pi(vd),Phi(va,:),Phi(vd,:),Sw(va),Sw(vd),Cp,CMP,sQm2,sQd2,QQ]=Sat_fast2([Phi(va,:);Phi(vd,:);Phi(vb,:)],[Pi(va,:);Pi(vd,:);Pi(vb,:)],Cp,[Sw(va,:);Sw(vd,:)],[Sw0(va,:);Sw0(vd,:)],RC,KWOG,KWOG_GY,WELL,BFRACM,...
    %       nw,Qz(:,ft+1),dt,Qf,GEOM,CMP,BB,TRM,M2FR,PR,fp,BXYZ,QQBND,QQoBND,CR_rc,sQm1,sQd1,QQ,ft,GYData);
       [Pi,Sw,Pw(:,ft),TL,DL,Phi,CMP,sQm2,sQc2,sQg2,sQd2,QQ,QQBND,QQoBND]=PressureCalcF2(Pi,Sw,Phi,Sw0,Pw(:,ft),Cp,TRM,KWOG,KWOG_GY,CMP,RC,WELL,fp,VEC,GEOM,DATA,M2FR,ft,PR,BXYZ,dt,Qf,GYData,BB,QQ,QQBND,QQoBND,ndt);
    else
       [Pi,Sw,Pw(:,ft),TL,DL,Phi,CMP,sQm2,sQd2,sQc2,sQg2,QQ,QQBND,QQoBND]=PressureCalc2(Pi,Sw,Phi,Sw0,Pwt(:,t),Cp,TRM,KWOG,KWOG_GY,CMP,RC,WELL,fp,VEC,GEOM,DATA,M2FR,ft,PR,BXYZ,dt,Qf,GYData,BB,QQ,QQBND,QQoBND); 
    end;
    So = 1 - Sw;
    MSw(:,1)=Sw(VEC.va);
    CSw(:,t+1)=Sw(VEC.vc);
    GSw(:,t+1)=Sw(VEC.vg);
    DSw(:,1)=Sw(VEC.vd);
    
    MSo(:,1)=So(VEC.va);
    CSo(:,t+1)=So(VEC.vc);
    GSo(:,t+1)=So(VEC.vg);
    DSo(:,1)=So(VEC.vd);
    
    MCp(:,1)=Cp(VEC.va);
    CCp(:,t+1)=Cp(VEC.vc);
    GCp(:,t+1)=Cp(VEC.vg);
    DCp(:,1)=Cp(VEC.vd);

    Qm(:,:,t+1) = sQm2;
%    Qm(CR_rc(1,1).won(:,3),:,t+1)=sQm1(CR_rc(1,1).won(:,3),:);
    Qc(:,:,t+1) = sQc2;
    Qg(:,:,t+1) = sQg2;   
    Qd(:,:,t+1) = sQd2;
%    Qd(CR_rc(1,2).won(:,3),:,t+1)=sQd1(CR_rc(1,2).won(:,3),:);   

    PpW(WELL.Won(:,3),t+1)=Pi(WELL.Won(:,1));
    Pwt(:,t+1)=Pw(:,ft);
    Ti(:,1)=Ti(:,1);
    
    if mod(t,zapt)==0
        if Ta/dt<=100
            zapt=1;
        elseif Ta/dt<=1000
            zapt=10;
        elseif Ta/dt<=1e4
            zapt=100;
        else
            zapt=1000;
        end;
        
        j=j+1;
        Pj(:,j)=Pi;
        Swj(:,j)=MSw;
        DSwj(:,j)=DSw;
        MCpj(:,j)=MCp;
        Tj(:,j)=Ti;
        waitbar(st/Ta)
    end;
%    sQo=sum(Qm(:,2,t+1)-Qm(:,3,t+1)+Qm(:,1,t+1))+sum(Qc(:,2,t+1)-Qc(:,3,t+1)+Qc(:,1,t+1))...
 %       +sum(Qg(:,2,t+1)-Qg(:,3,t+1)+Qg(:,1,t+1))+sum(Qd(:,2,t+1)-Qd(:,3,t+1)+Qd(:,1,t+1));
    %sum(sQo(:))
    % dQ(t)=sum(Sw0.*[dV;dVC;dVG])-sum([Sw;Cw(:,t+1);Gw(:,t+1)].*[dV;dVC;dVG])-sum(sQo(:));
    
    %%prob.progress;
 %   c_lik=1-Qm(:,3,t+1)./Qm(:,2,t+1);
 %   c_lik(isnan(c_lik)==1)=0;
 %   Uf(:,ft+1:end)=Uf(:,ft+1:end).*repmat(c_lik<PwQC_bnd(:,7),1,size(Uf(:,ft+1:end),2));
 %   qin=-Qm(:,1,t+1)/dt;
 %   qo=Qm(:,3,t+1)/dt;
 %   qo(Uf(:,ft+1)==-1)=inf;
 %   Qmin=repmat(PwQC_bnd(:,8),1,1);
 %   Uf(:,ft+1:end)=Uf(:,ft+1:end).*repmat(qo>Qmin,1,size(Uf(:,ft+1:end),2));
    
    dt1(t+1)=dt;
    st=st+dt;
    t_flag=st<Ta;
end;

j=j+1;
Pj(:,j)=Pi;
Swj(:,j)=MSw;
DSwj(:,j)=DSw;
CSwj(:,j)=Sw(vc);
Soj(:,j)=MSo;
DSoj(:,j)=DSo;
MCpj(:,j)=Cp(va);
Tj(:,j)=Ti;
waitbar(st/Ta)
Q = 0;
Uf = WELL.Uf;
NDT = 1;
dtz = 1;
[Q,Pw,PpW,dtz]=Q2Sut(Qm,Qc,Qg,Qd,Pwt(:,1:t+1),PpW(:,1:t+1),dt1,Ta);
% GY_Pxy(p)=GY_Pxy;
% save('GY_Pxy.mat','GY_Pxy')
%  Bnd_xy(p)=DATA.BndXY;
%  GY_Pxy(Bnd_xy==1)'
toc

