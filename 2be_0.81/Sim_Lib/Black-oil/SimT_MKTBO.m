function [Pj,Pbj,Swj,Soj,Tj,MCpj,p,Q,Pw,PpW,CSwj,NDT,Uf,dt1,dV0,ka,dtz,DSwj]=SimT_MKTBO(PR,C,A2C,GData,BB,A2B,DData,dVC,dVB,DATA,WData,GYData,fll,CR_GRUP)
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

[G,M2FR.A2G,dVG,WELL.WonG,GEOM.Lg,Mp_g]=Rasp(GData,3);   %Распаковка структуры для гор. трещ.
[D,M2FR.A2D,dVD,WELL.WonD,GEOM.Ld,Mp_d]=Rasp(DData,4);   %Распаковка структуры для двойной среды

WELL.Won=DATA.Won;       %номера ячеек скважин, 
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

M2FR.A2C(raz+1:end+raz,:)=M2FR.A2C;  %A2C - связь пор с вертикальными трещинами
M2FR.A2C(1:raz,:)=0;
M2FR.A2C=M2FR.A2C(p,:);
M2FR.A2G=M2FR.A2G(p,:);  % - с горизонтальными
M2FR.A2D=M2FR.A2D(p,:);  % - вторая среда (мелкие трещины)
M2FR.A2B=M2FR.A2B(p,:);

for i=1:size(WELL.Won,1)
    WELL.Won(i,1)=find(WELL.Won(i,1)==p);
end;

na=size(XY,1);  nc=size(C,1);  ng=size(G,1);  nd=size(D,1); nb=size(BB,1); nw=size(Qz,1);
NCELL.na=na;  NCELL.nc=nc;  NCELL.ng=ng;  NCELL.nd=nd; NCELL.nb=nb; NCELL.nw=nw;
% na - число узлов в поровой матрице, nc - число узлов в трещинах, ng - число узлов в гориз. трещинах, 
% nd - число узлов в двойной среде, nb - число узлов в аквифере

va=1:na;                                VEC.va=va; % - поровая матрица
vc=na+1:na+nc;                          VEC.vc=vc; % - трещины
vg=na+nc+1:na+nc+ng;                    VEC.vg=vg; % - горизонтальные трещины
vd=na+nc+ng+1:na+nc+ng+nd;              VEC.vd=vd; % - двойная среда
vb=na+nc+ng+nd+1:na+nc+ng+nd+nb;        VEC.vb=vb; % - аквифер
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
    GEOM.Ke,GEOM.L,B,S,H1,GEOM.K(:,1),DData.Kd,Ke_gy,BndXY(p),BndZ(p),nb,dVCG);  %формирование геометрической проводимостей
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

Pb = 70*ones(na+nc+ng+nd,1);

Ti(va,1)=T;
Ti(vc,1)=T(RC.ACr);
Ti(vg,1)=T(RC.AGr);
Ti(vd,1)=T(RC.ADr);
%Ti(vb,1)=GYData.T0;

Phi(VEC.va,:) = [Pi(VEC.va) - Z(VEC.va)*Ro(1)*PR.g, Pi(VEC.va) - Z(VEC.va)*Ro(2)*PR.g, (Pi(VEC.va) - Z(VEC.va)*Ro(3)*PR.g)];
Phi(VEC.vc,:) = [Pi(VEC.vc) - Z(RC.ACr)*Ro(1)*PR.g, Pi(VEC.vc) - Z(RC.ACr)*Ro(2)*PR.g, (Pi(VEC.vc) - Z(RC.ACr)*Ro(3)*PR.g)];
Phi(VEC.vg,:) = [Pi(VEC.vg) - Z(RC.AGr)*Ro(1)*PR.g, Pi(VEC.vg) - Z(RC.AGr)*Ro(2)*PR.g, (Pi(VEC.vg) - Z(RC.AGr)*Ro(3)*PR.g)];
Phi(VEC.vd,:) = [Pi(VEC.vd) - Z(RC.ADr)*Ro(1)*PR.g, Pi(VEC.vd) - Z(RC.ADr)*Ro(2)*PR.g, (Pi(VEC.vd) - Z(RC.ADr)*Ro(3)*PR.g)];
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

CMP.B0 = Bo;
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
 
 VSAT.vg = VEC.v(Pb(VEC.v) >= Pi(VEC.v)); VSAT.vp = VEC.v(Pb(VEC.v) < Pi(VEC.v)); 
 
dt = 0.25; 
Pwt(:,1)=Pw(:,1);
 
[CMP]=SGimBO0(Sw,So,zc,PR.rs,Pi,Pb,CMP,VSAT,dt,GEOM.dV);

TL = ones(na,na);
DL = ones(nd,nd);

 BXYZ.mxy = find(TRM.Txyz_GY_A(:,1));
 BXYZ.dxy = find(TRM.Txyz_GY_D(:,1));
 BXYZ.mz = find(TRM.Txyz_GY_A(:,2));
 BXYZ.dz = find(TRM.Txyz_GY_A(:,2)); 
       
while t_flag==1
    t=t+1;  
%    if t==1 && dtt==0
%        dt=1;
%        dtt=1;
%    else
%        dtt=dtt1;
%    end;
    
    ft=floor(st)+1;
    % if ~isempty(find(CpW(:,t)~=0)) fp=1; else fp=0; end;
    fp=0;
    Qf=Qz(:,ft);
    
    Sw0 = Sw;
    So0 = So;
    CMP.Rs(:,1) = CMP.Rs(:,2);
    CMP.Bw(:,1) = CMP.Bw(:,2);
    CMP.Bo(:,1) = CMP.Bo(:,2);
    CMP.Bg(:,1) = CMP.Bg(:,2);
    CMP.Mp(:,1) = CMP.Mp(:,2);
    
    KWOG.w(1:na,1)=Sat_cal(MSw,MSo,1,1,as,aw); %water
    KWOG.o(1:na,1)=Sat_cal(MSw,MSo,2,1,as,aw); %oil 
    KWOG.g(1:na,1)=Sat_cal(MSw,MSo,3,1,as,aw); %gas  ///////////
    
    KWOG.w(na+1:na+nc+ng,1)=Sat_cal([CSw(:,t);GSw(:,t)],[CSo(:,t);GSo(:,t)],1,1,ts,tw); %water
    KWOG.o(na+1:na+nc+ng,1)=Sat_cal([CSw(:,t);GSw(:,t)],[CSo(:,t);GSo(:,t)],2,1,ts,tw); %oil
    KWOG.g(na+1:na+nc+ng,1)=Sat_cal([CSw(:,t);GSw(:,t)],[CSo(:,t);GSo(:,t)],3,1,ts,tw); %gas  ////////////

   [W1,~,~,~,~,QQ.QQm,QQ.QQmwog]=Well_MKTBO0(Pi(VEC.va),Pwt(WELL.Won(:,3),t),WELL.Won,WELL.Uf(WELL.Won(:,3)),Cp(VEC.va),PR,WELL.CpW(WELL.Won(:,3)),CMP,KWOG,VEC.va,na);% W1 - проводимость по всей жидкости, W6 - только для воды, W7 - полимер 
   [W1C,~,~,~,~,QQ.QQc,QQ.QQcwog]=Well_MKTBO0(Pi(VEC.vc),Pwt(WELL.WonC(:,3),t),WELL.WonC,WELL.Uf(WELL.WonC(:,3)),Cp(VEC.vc),PR,WELL.CpW(WELL.WonC(:,3)),CMP,KWOG,VEC.vc,nc);
   [W1G,~,~,~,~,QQ.QQg,QQ.QQgwog]=Well_MKTBO0(Pi(VEC.vg),Pwt(WELL.WonG(:,3),t),WELL.WonG,WELL.Uf(WELL.WonG(:,3)),Cp(VEC.vg),PR,WELL.CpW(WELL.WonG(:,3)),CMP,KWOG,VEC.vg,ng);
   [W1D,~,~,~,~,QQ.QQd,QQ.QQdwog]=Well_MKTBO0(Pi(VEC.vd),Pwt(WELL.WonD(:,3),t),WELL.WonD,WELL.Uf(WELL.WonD(:,3)),Cp(VEC.vd),PR,WELL.CpW(WELL.WonD(:,3)),CMP,KWOG,VEC.vd,nd);
   
   [WBND,QQBND,QQwoBND,KWOG_GY]=GY_bildBO0(GYData,Pi([VEC.va,VEC.vd],1),Sw([VEC.va,VEC.vd],1),So([VEC.va,VEC.vd],1),Cp([VEC.va,VEC.vd],1),RC,TRM.Txyz_GY_A,TRM.Txyz_GY_D,PR,CMP,BXYZ,VEC);
  
   %dt=vibor_t2(dtt,Pi,RC,dVCG,TL,W1,WELL,Pwt(:,t),na,PR,st,Ta,Sw,Sw0,dt,Nl,va,vd,DL,W1D,nd,dV1,dV2);   
 
   [Pi,Pb,Sw,So,Pw(:,ft),TL,DL,Phi,CMP,sQm2,sQc2,sQg2,sQd2,VSAT,QQ,QQBND,QQwoBND]=PressureCalcBO(Pi,Sw,So,Phi,Sw0,So0,Pb,Pwt(:,t),Cp,TRM,KWOG,KWOG_GY,CMP,RC,WELL,fp,VEC,GEOM,DATA,NCELL,M2FR,ft,PR,BXYZ,dt,Qf,VSAT,GYData,BB,QQ,QQBND,QQwoBND);
 
 %   if isempty(RC.Cr)==0 || isempty(RC.Gr)==0
 %       Sw = Sw0;
 %       So = So0;
 %       [ndt,~,~]=vibor_t(ndt,0,Pi(VEC.vc),Pi(VEC.vg),Pw(:,ft+1),PR,RC,dt,Sw0([VEC.vc,VEC.vg]),Sw([VEC.vc,VEC.vg]),...
 %           CL(RC.Cr2+(RC.Cc2-1)*nc),GL(RC.Gr2+(RC.Gc2-1)*ng),dVCG(VEC.vc),dVCG(VEC.vg),W1C,W1G,WELL.WonC,WELL.WonG,1,0);
        
  %      [BFRACM,Sw(VEC.vc),Sw(VEC.vg),So(VEC.vc),So(VEC.vg),ndt,sQ1,sQ2,sQm1,sQd1,QWOG]=fun3(RC,Pi,Phi,Sw,So,Cp,PR,TRM,M2FR,WELL,Pw,dt,CR_rc,Qz(:,ft+1),Qf,ndt,Pi0,GEOM,DATA.Lc,P0,SGM,ft,QWOG,VSAT);
        
   %     [Sw,So,Cp,sQm2,sQd2,QWOG]=Sat_fast_3([Phi(VEC.va,:);Phi(VEC.vd,:);Phi(VEC.vb,:)],Pi,Cp,Sw,So,Sw0,So0,RC,KWOG,WELL,BFRACM,...
   %         nw,Qz(:,ft+1),dt,Qf,GEOM,SGMMD,BB,TRM,M2FR,PR,fp,BXYZ,QBOUND,CR_rc,sQm1,sQd1,QWOG,ft,VSAT,GYData);
%    end;
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
    Qc(:,:,t+1) = sQc2;
    Qg(:,:,t+1) = sQg2;    
    Qd(:,:,t+1) = sQd2;
 
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
        Pbj(:,j)=Pb;
        Swj(:,j)=MSw;
        DSwj(:,j)=DSw;
        MCpj(:,j)=MCp;
        Tj(:,j)=Ti;
        waitbar(st/Ta)
    end;
    sQo=sum(Qm(:,2,t+1)-Qm(:,3,t+1)+Qm(:,1,t+1))+sum(Qc(:,2,t+1)-Qc(:,3,t+1)+Qc(:,1,t+1))...
        +sum(Qg(:,2,t+1)-Qg(:,3,t+1)+Qg(:,1,t+1))+sum(Qd(:,2,t+1)-Qd(:,3,t+1)+Qd(:,1,t+1));
    %sum(sQo(:))
    % dQ(t)=sum(Sw0.*[dV;dVC;dVG])-sum([Sw;Cw(:,t+1);Gw(:,t+1)].*[dV;dVC;dVG])-sum(sQo(:));

    dt1(t+1)=dt;
    st=st+dt;
    t_flag=st<Ta;
end;

j=j+1;
Pj(:,j)=Pi;
Pbj(:,j)=Pb;
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

