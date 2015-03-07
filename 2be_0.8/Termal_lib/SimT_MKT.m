function [Pj,Swj,Soj,Tj,MCpj,p,Q,Pw,PpW,CSwj,NDT,Uf,dt1,dV0,ka,dtz,DSwj]=SimT_MKT(PR,C,A2C,GData,BB,A2B,DData,dVC,dVB,DATA,WData,GYData,fll,CR_GRUP)
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
Pi(vb,1)=GYData.P0;

Ti(va,1)=T;
Ti(vc,1)=T(RC.ACr);
Ti(vg,1)=T(RC.AGr);
Ti(vd,1)=T(RC.ADr);
Ti(vb,1)=GYData.T0;


Phi(VEC.va,:) = [Pi(VEC.va) - Z(VEC.va)*Ro(1)*PR.g, Pi(VEC.va) - Z(VEC.va)*Ro(2)*PR.g, 0*(Pi(VEC.va) - Z(VEC.va)*Ro(3)*PR.g)];
Phi(VEC.vc,:) = [Pi(VEC.vc) - Z(RC.ACr)*Ro(1)*PR.g, Pi(VEC.vc) - Z(RC.ACr)*Ro(2)*PR.g, 0*(Pi(VEC.vc) - Z(RC.ACr)*Ro(3)*PR.g)];
Phi(VEC.vg,:) = [Pi(VEC.vg) - Z(RC.AGr)*Ro(1)*PR.g, Pi(VEC.vg) - Z(RC.AGr)*Ro(2)*PR.g, 0*(Pi(VEC.vg) - Z(RC.AGr)*Ro(3)*PR.g)];
Phi(VEC.vd,:) = [Pi(VEC.vd) - Z(RC.ADr)*Ro(1)*PR.g, Pi(VEC.vd) - Z(RC.ADr)*Ro(2)*PR.g, 0*(Pi(VEC.vd) - Z(RC.ADr)*Ro(3)*PR.g)];

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

SGM.Bwog(:,1:2)=Bo(1)*ones(na+nc+ng+nd+nb,2);       %/////////
SGM.Bwog(:,3:4)=Bo(2)*ones(na+nc+ng+nd+nb,2);       %/////////
SGM.Bwog(:,5:6)=Bo(3)*ones(na+nc+ng+nd+nb,2);
SGM.Mp=repmat([Mp.*ones(na,1);Mp_c;Mp_g;Mp_d;Mp_b],1,2);
GEOM.dV = dVCG./SGM.Mp(:,1); 

Sw=[MSw(:,1);CSw(:,1);GSw(:,1);DSw(:,1);BSw(:,1)];
So=[MSo(:,1);CSo(:,1);GSo(:,1);DSo(:,1);BSo(:,1)];  %//////////
Cp=[MCp(:,1);CCp(:,1);GCp(:,1);DCp(:,1);BCp(:,1)];

V0=sum(dVCG.*(1-Sw));
P0=Pi;

dtt1=dtt;
Wf0=WELL.Won(:,2);
zapt=1;

dt = 1; 
[SGM]=SGim0(Sw,So,zc,PR.rs,Pi,1,P0,SGM,GEOM,VEC,dt);

[W1,W6,Wo,Wg,W7]=Well_MKT(WELL.Won,WELL.Uf(WELL.Won(:,3),1),Sw(VEC.va,1),Cp(VEC.va,1),PR.aw,PR.as,PR,WELL.CpW(WELL.Won(:,3),1),SGM);% W1 - проводимость по всей жидкости, W6 - только для воды, W7 - полимер 
[W1C,W6C,WoC,WgC,W7C]=Well_MKT(WELL.WonC,WELL.Uf(WELL.WonC(:,3),1),Sw(VEC.vc,1),Cp(VEC.vc,1),PR.tw,PR.ts,PR,WELL.CpW(WELL.WonC(:,3),1),SGM);
[W1G,W6G,WoG,WgG,W7G]=Well_MKT(WELL.WonG,WELL.Uf(WELL.WonG(:,3),1),Sw(VEC.vg,1),Cp(VEC.vg,1),PR.tw,PR.ts,PR,WELL.CpW(WELL.WonG(:,3),1),SGM);
[W1D,W6D,WoD,WgD,W7D]=Well_MKT(WELL.WonD,WELL.Uf(WELL.WonD(:,3),1),Sw(VEC.vd,1),Cp(VEC.vd,1),PR.tw,PR.ts,PR,WELL.CpW(WELL.WonD(:,3),1),SGM);

[QWOG.QQm,QQmwo]=QIter(zeros(size(WELL.Won,1),1),W6,Wo,Wg,Pi(VEC.va),WELL.Won(:,1),Pw(WELL.Won(:,3),1),SGM);
[QWOG.QQc,QQcwo]=QIter(zeros(size(WELL.WonC,1),1),W6C,WoC,WgC,Pi(VEC.vc),WELL.WonC(:,1),Pw(WELL.WonC(:,3),1),SGM);
[QWOG.QQg,QQgwo]=QIter(zeros(size(WELL.WonG,1),1),W6G,WoG,WgG,Pi(VEC.vg),WELL.WonG(:,1),Pw(WELL.WonG(:,3),1),SGM);
[QWOG.QQd,QQdwo]=QIter(zeros(size(WELL.WonD,1),1),W6D,WoD,WgD,Pi(VEC.vd),WELL.WonD(:,1),Pw(WELL.WonD(:,3),1),SGM);     

QWOG.QQmwo(:,1) = sparse(WELL.Won(:,1),ones(1,size(WELL.Won,1)),QQmwo(:,1),na,1);
QWOG.QQcwo(:,1) = sparse(WELL.WonC(:,1),ones(1,size(WELL.WonC,1)),QQcwo(:,1),nc,1);
QWOG.QQgwo(:,1) = sparse(WELL.WonG(:,1),ones(1,size(WELL.WonG,1)),QQgwo(:,1),ng,1);
QWOG.QQdwo(:,1) = sparse(WELL.WonD(:,1),ones(1,size(WELL.WonD,1)),QQdwo(:,1),nd,1);
   
QWOG.QQmwo(:,2) = sparse(WELL.Won(:,1),ones(1,size(WELL.Won,1)),QQmwo(:,2),na,1);
QWOG.QQcwo(:,2) = sparse(WELL.WonC(:,1),ones(1,size(WELL.WonC,1)),QQcwo(:,2),nc,1);
QWOG.QQgwo(:,2) = sparse(WELL.WonG(:,1),ones(1,size(WELL.WonG,1)),QQgwo(:,2),ng,1);
QWOG.QQdwo(:,2) = sparse(WELL.WonD(:,1),ones(1,size(WELL.WonD,1)),QQdwo(:,2),nd,1);

while t_flag==1
    t=t+1;

    if t==1 && dtt==0
        dt=1;
        dtt=1;
    else
        dtt=dtt1;
    end;
    
    Sw0 = Sw;
    So0 = So;
    Pi0=Pi;
    SGM.Rs(:,1) = SGM.Rs(:,2);
    SGM.Bwog(:,1) = SGM.Bwog(:,2);
    SGM.Bwog(:,3) = SGM.Bwog(:,4);
    SGM.Bwog(:,5) = SGM.Bwog(:,6);
    SGM.Mp(:,1) = SGM.Mp(:,2);

    ft=floor(st);
    % if ~isempty(find(CpW(:,t)~=0)) fp=1; else fp=0; end;
    fp=0;
    Qf=Qz(:,ft+1);
  
    KWOG.w(1:na,1)=Sat_cal(MSw,1,1,as,aw); %water
    KWOG.o(1:na,1)=Sat_cal(MSw,2,1,as,aw); %oil 
    KWOG.g(1:na,1)=Sat_cal(MSw,3,1,as,aw); %gas  ///////////
    
    KWOG.w(na+1:na+nc+ng,1)=Sat_cal([CSw(:,t);GSw(:,t)],1,1,ts,tw); %water
    KWOG.o(na+1:na+nc+ng,1)=Sat_cal([CSw(:,t);GSw(:,t)],2,1,ts,tw); %oil
    KWOG.g(na+1:na+nc+ng,1)=Sat_cal([CSw(:,t);GSw(:,t)],3,1,ts,tw); %gas  ////////////
    
    [BOUNDFL]=GY_bild(GYData,Pi([VEC.va,VEC.vd],1),Sw([VEC.va,VEC.vd],1),Cp([VEC.va,VEC.vd],1),RC,TRM.Txyz_GY_A,TRM.Txyz_GY_D,PR);

    QBOUND.Qm = BOUNDFL.b1gm(:,2).*(GYData.GY_Pz - Pi(VEC.va))+ BOUNDFL.b1gm(:,1).*(GYData.GY_Pxy - Pi(VEC.va));
    QBOUND.Qc = BOUNDFL.b1gc;
    QBOUND.Qg = BOUNDFL.b1gg;
    QBOUND.Qd = BOUNDFL.b1gd(:,2).*(GYData.GY_Pz(RC.ADr) - Pi(VEC.vd))+ BOUNDFL.b1gd(:,1).*(GYData.GY_Pxy(RC.ADr) - Pi(VEC.vd));
    QBOUND.Qb = BOUNDFL.b1gb.*(GYData.P0 - Pi(VEC.vb))+ BOUNDFL.b1gb(:,1).*(GYData.P0 - Pi(VEC.vb));

    [Pi,CL,GL,QWOG,Phi,SGMMD,sQm2,sQd2]=PressureCalc(Pi,Sw,So,Phi,Sw0,So0,Cp,TRM,KWOG,SGM,RC,WELL,fp,VEC,GEOM,DATA,NCELL,M2FR,ft,PR,QWOG,BOUNDFL,dt,t,Qf,QBOUND,Qz);
    
    if isempty(RC.Cr)~=0 || isempty(RC.Gr)~=0
        [ndt,~,~]=vibor_t(ndt,0,Pi(VEC.vc),Pi(VEC.vg),Pw(:,ft+1),PR,RC,dt,Sw0([VEC.vc,VEC.vg]),Sw([VEC.vc,VEC.vg]),...
            CL(RC.Cr2+(RC.Cc2-1)*nc),GL(RC.Gr2+(RC.Gc2-1)*ng),dVCG(VEC.vc),dVCG(VEC.vg),W1C,W1G,WELL.WonC,WELL.WonG,1,0);
        
        [BFRACM,Sw(VEC.vc),Sw(VEC.vg),So(VEC.vc),So(VEC.vg),ndt,sQ1,sQ2,sQm1,sQd1,QWOG]=fun3(RC,Pi,Phi,Sw,So,Cp,PR,TRM,M2FR,WELL,Pw,dt,CR_rc,Qz(:,ft+1),Qf,ndt,Pi0,GEOM,DATA.Lc,P0,SGM,ft,QWOG);
        
        [Sw,So,Cp,sQm2,sQd2,QWOG]=Sat_fast_3([Phi(VEC.va,:);Phi(VEC.vd,:);Phi(VEC.vb,:)],Pi,Cp,Sw,So,Sw0,So0,RC,KWOG,WELL,BFRACM,...
            nw,Qz(:,ft+1),dt,Qf,GEOM,SGMMD,BB,TRM,M2FR,PR,fp,BOUNDFL,QBOUND,CR_rc,sQm1,sQd1,QWOG,ft);
    end;
    MSw(:,1)=Sw(VEC.va);
    CSw(:,t+1)=Sw(VEC.vc);
    GSw(:,t+1)=Sw(VEC.vg);
    DSw(:,1)=Sw(VEC.vd);
    
    MCp(:,1)=Cp(VEC.va);
    CCp(:,t+1)=Cp(VEC.vc);
    GCp(:,t+1)=Cp(VEC.vg);
    DCp(:,1)=Cp(VEC.vd);

    Qm(:,:,t+1) = sQm2;
    Qm(CR_rc(1,1).won(:,3),:,t+1)=sQm1(CR_rc(1,1).won(:,3),:);
    Qc(:,:,t+1) = sQ1;
    Qg(:,:,t+1) = sQ2;   
    Qd(:,:,t+1) = sQd2;
    Qd(CR_rc(1,2).won(:,3),:,t+1)=sQd1(CR_rc(1,2).won(:,3),:);   

    PpW(WELL.Won(:,3),t+1)=Pi(WELL.Won(:,1));
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
    sQo=sum(Qm(:,2,t+1)-Qm(:,3,t+1)+Qm(:,1,t+1))+sum(Qc(:,2,t+1)-Qc(:,3,t+1)+Qc(:,1,t+1))...
        +sum(Qg(:,2,t+1)-Qg(:,3,t+1)+Qg(:,1,t+1))+sum(Qd(:,2,t+1)-Qd(:,3,t+1)+Qd(:,1,t+1));
    %sum(sQo(:))
    % dQ(t)=sum(Sw0.*[dV;dVC;dVG])-sum([Sw;Cw(:,t+1);Gw(:,t+1)].*[dV;dVC;dVG])-sum(sQo(:));
    
    %%prob.progress;
    c_lik=1-Qm(:,3,t+1)./Qm(:,2,t+1);
    c_lik(isnan(c_lik)==1)=0;
    Uf(:,ft+1:end)=Uf(:,ft+1:end).*repmat(c_lik<PwQC_bnd(:,7),1,size(Uf(:,ft+1:end),2));
    qin=-Qm(:,1,t+1)/dt;
    qo=Qm(:,3,t+1)/dt;
    qo(Uf(:,ft+1)==-1)=inf;
    Qmin=repmat(PwQC_bnd(:,8),1,1);
    Uf(:,ft+1:end)=Uf(:,ft+1:end).*repmat(qo>Qmin,1,size(Uf(:,ft+1:end),2));
    
    st=st+dt;
    dt1(t)=dt;
    dt=vibor_t2(dtt,Pi,RC,dVCG,TL,W1,Won,Pw(:,ft+1),na,PR,st,Ta,Sw,Sw0,dt,Nl,WonM,va,vd,DL,W1D,WonD,nd,dV1,dV2);
    t_flag=st~=Ta;
    
end;

j=j+1;
Pj(:,j)=Pi;
Swj(:,j)=MSw;
DSwj(:,j)=DSw;
MCpj(:,j)=MCp;
Tj(:,j)=Ti;
waitbar(st/Ta)
[Q,Pw,PpW,dtz]=Q2Sut(Qm,Qc,Qg,Qd,Pwt(:,1:t+1),PpW(:,1:t+1),dt1,Ta);
% GY_Pxy(p)=GY_Pxy;
% save('GY_Pxy.mat','GY_Pxy')
%  Bnd_xy(p)=DATA.BndXY;
%  GY_Pxy(Bnd_xy==1)'
toc

% 
% 
%     kj=0;
%     flag_gim=1;
%     flag_pwq=1;
%     
%     while flag_gim==1 && kj<10
%         kj=kj+1;
%         
%         [Clp,Cwp,Cws,A,Bwo,Mp]=SGim(dVCG./Mp(:,1),Sw,Mp,zc,Bwo,Pi,1,P0,va,vc,vg,vd,vb,dt);
%         
%         [TL,TW,TP]=Potok_MKT(TTM,Pi(va,1),kfw(va),kfo(va),MCp(:,1),mu,RC.Arc2,mup,fp,kms(1),L,Ke,Ro,A(va),dZ(1,:));
%         [CL,CW,~]=Potok_Tube(TC,Pi(vc,1),CSw(:,t),CCp(:,t),PR,mup,fp,kms(2),DATA.Lc,RC.Cr2,RC.Cc2,nc,A(vc),dZ(2,:));
%         [GL,GW,~]=Potok_Tube(TG,Pi(vg,1),GSw(:,t),GCp(:,t),PR,mup,fp,kms(3),Lg,RC.Gr2,RC.Gc2,ng,A(vg),dZ(3,:));
%         [DL,DW,DP]=Potok_Tube(TD,Pi(vd,1),DSw(:,1),DCp(:,1),PR,mup,fp,kms(4),Ld,RC.Dr2,RC.Dc2,nd,A(vd),dZ(4,:));
%         
%          %qe=ciklik(TL,50);
%          
%         [Gr,Grw]=Gravity(TL,TW,CL,CW,GL,GW,DL,DW,[],A,dZ);
%         
%         [A2CL,~,~]=Obmen_T2M(A2C,Pi(va,1),Pi(vc,1),MSw(:,1),CSw(:,t),K(:,1),PR,MCp(:,1),CCp(:,t));
%         [A2GL,~,~]=Obmen_T2M(A2G,Pi(va,1),Pi(vg,1),MSw(:,1),GSw(:,t),K(:,1),PR,MCp(:,1),GCp(:,t));
%         
%         [A2DL,A2DW,A2DP]=Obmen_T2M(A2D,Pi(va,1),Pi(vd,1),MSw(:,1),DSw(:,1),K(:,1),PR,MCp(:,1),DCp(:,1));
%         [A2BL,A2BW,A2BP]=Obmen_T2M(A2B,Pi(va,1),Pi(vb,1),MSw(:,1),BSw(:,1),ones(na,1),PR,MCp(:,1),BCp(:,1));
%         [D2BL,D2BW,D2BP]=Obmen_T2M(D2B,Pi(vd,1),Pi(vb,1),DSw(:,1),BSw(:,1),K(:,1),PR,DCp(:,1),BCp(:,1));
%         
%         [D2CL,~,~]=Obmen_T2M(C2D,Pi(vd,1),Pi(vc,1),DSw(:,1),CSw(:,t),K(:,1),PR,DCp(:,1),CCp(:,t));
%         [D2GL,~,~]=Obmen_T2M(G2D,Pi(vd,1),Pi(vg,1),DSw(:,1),GSw(:,t),K(:,1),PR,DCp(:,1),GCp(:,t));
%         
%         %Wf=Wf0.*(1-0.0001*(P0(Won)-Pi(Won))).^3;
%         [W1,W6,W7]=Well_MKT(Wf,Won(:,1),Uf(WonM,ft+1),MSw(:,1),MCp(:,1),aw,as,mu,mup,CpW(WonM,ft+1),A(va));
%         [W1C,W6C,W7C]=Well_MKT(WonC(:,2),WonC(:,1),Uf(WonC(:,3),ft+1),CSw(:,t),CCp(:,t),tw,ts,mu,mup,CpW(WonC(:,3),ft+1),A(vc));
%         [W1G,W6G,W7G]=Well_MKT(WonG(:,2),WonG(:,1),Uf(WNG,ft+1),GSw(:,t),GCp(:,t),tw,ts,mu,mup,CpW(WNG,ft+1),A(vg));
%         [W1D,W6D,W7D]=Well_MKT(WonD(:,2),WonD(:,1),Uf(WonD(:,3),ft+1),DSw(:,1),DCp(:,1),tw,ts,mu,mup,CpW(WonD(:,3),ft+1),A(vd));
%         
%         A1=TL-sparse(Won(:,1),Won(:,1),W1,na,na)-sparse(1:na,1:na,sum(A2CL,2)+sum(A2GL,2)+sum(A2BL,2)+sum(A2DL,2)+Clp(va)+sum(b1gm(:,1:2),2),na,na);  %Матрица коэф. для пор
%         C1=CL-sparse(1:nc,1:nc,sum(A2CL,1)+sum(D2CL,1)+Clp(vc)',nc,nc)-sparse(WonC(:,1),WonC(:,1),W1C,nc,nc);                       %Матрица коэф. для вертикальных трещ.
%         G1=GL-sparse(1:ng,1:ng,sum(A2GL,1)+sum(D2GL,1)+Clp(vg)',ng,ng)-sparse(WonG(:,1),WonG(:,1),W1G,ng,ng);                       %Матрица коэф. для гориз. трещ.
%         D1=DL-sparse(WonD(:,1),WonD(:,1),W1D,nd,nd)-sparse(1:nd,1:nd,sum(A2DL,1)'+sum(D2BL,2)+sum(D2CL,2)+sum(D2GL,2)+Clp(vd)+sum(b1gd(:,1:2),2),nd,nd);                       %Матрица коэф. для двойной пор.
%         B1=BB-sparse(1:nb,1:nb,sum(A2BL,1)+sum(D2BL,1)+Clp(vb)'+b1gb',nb,nb);                                                             %Матрица коэф. для границ
%         
%         W2M=sparse(WonM,Won(:,1),W1,nw,na);
%         W2C=sparse(WonC(:,3),WonC(:,1),W1C,nw,nc);
%         W2G=sparse(WNG,WonG(:,1),W1G,nw,ng);
%         W2D=sparse(WonD(:,3),WonD(:,1),W1D,nw,nd);
%         W2B=sparse(nw,nb);
% 
%         WM1=[W2M,W2C,W2G,W2D,W2B];
%         WM2=WM1';
%         W3vec=sparse(WonM,1,W1,nw,1)+sparse(WonC(:,3),1,W1C,nw,1)+sparse(WNG,1,W1G,nw,1)+sparse(WonD(:,3),1,W1D,nw,1);
%         WM3=-sparse(1:nw,1:nw,W3vec,nw,nw);
%         
%         AM=[A1,   A2CL, A2GL, A2DL, A2BL;
%             A2CL',  C1,  C2GL, D2CL', C2BL;
%             A2GL', C2GL', G1,  D2GL', G2BL;
%             A2DL', D2CL, D2GL,  D1,  D2BL;
%             A2BL', C2BL',G2BL',D2BL', B1];
%             
%         while  flag_pwq==1
%            
%             PwNl=repmat(Pw(:,ft+1),Nl,1);
%             % PwNl=PwNl(ka1==1);
%             
%             b1wm=sparse(Won(:,1),ones(1,size(Won,1)),-W1.*PwNl(WonM),na,1);
%             b1wc=sparse(WonC(:,1),ones(1,size(WonC,1)),-W1C.*Pw(WonC(:,3),ft+1),nc,1);
%             b1wg=sparse(WonG(:,1),ones(1,size(WonG,1)),-W1G.*Pw(WNG,ft+1),ng,1);
%             b1wd=sparse(WonD(:,1),ones(1,size(WonD,1)),-W1D.*PwNl(WonD(:,3)),nd,1);
%             
%             b1wm=b1wm.*(sum(W2M(Qf==0,:),1)~=0)';
%             b1wc=b1wc.*(sum(W2C(Qf==0,:),1)~=0)';
%             b1wg=b1wg.*(sum(W2G(Qf==0,:),1)~=0)';
%             b1wd=b1wd.*(sum(W2D(Qf==0,:),1)~=0)';
%             %  b1wm
% 
%             BM=[b1wm;b1wc;b1wg;b1wd;b1wb]+[-b1gm(:,2).*GY_Pz-b1gm(:,1).*GY_Pxy;b1gc;b1gg;...
%                 [-b1gd(:,2).*GY_Pz(vad)-b1gd(:,1).*GY_Pxy(vad)];-b1gb.*GY_Pxy2]-(Clp.*Pi)-Gr;
%             
%             BLGY_GIM=[[-b1gm(:,2).*GY_Pz-b1gm(:,1).*GY_Pxy];[-b1gd(:,2).*GY_Pz(vad)-b1gd(:,1).*GY_Pxy(vad)];...
%                 -b1gb.*GY_Pxy2]-(Clp([va,vd,vb]).*Pi([va,vd,vb]))-Gr([va,vd,vb]);
%             
%             W2M1=WM1(Qf~=0,:);
%             W2M2=WM2(:,Qf~=0);
%             W2M3=WM3(Qf~=0,Qf~=0);
% 
%             Pt=[BM',Qz(Qf~=0,ft+1)']/[AM,W2M2;W2M1,W2M3];
%            
%             flag_gim=sum(abs(Pt(1:na+nc+ng+nd+nb)-Pt0(1:na+nc+ng+nd+nb))./Pt(1:na+nc+ng+nd+nb)>=1e-6)~=0;
%             flag_gim=flag_gim*(sum(zc==0)==0);
%             
%             pw=Pw(:,ft+1);
%             pw(Qf~=0)=Pt(na+nc+ng+nd+nb+1:end);
%             
%             %[flag_pwq,Pw(:,ft+1),Qz(:,ft+1),Qf]=Chek_bond(pw,Pt(Won),W1,Uf(WonM,ft+1),Qf,PwQC_bnd);
%             
%             qm=QBild(W1,W6,W7,Pt(va)',Uf(WonM,ft+1),Won(:,1),dt,pw(WonM),WonM,nw);
%             qc=QBild(W1C,W6C,W7C,Pt(vc)',Uf(WonC(:,3),ft+1),WonC(:,1),dt,pw(WonC(:,3)),WonC(:,3),nw);
%             qg=QBild(W1G,W6G,W7G,Pt(vg)',Uf(WNG,ft+1),WonG(:,1),dt,pw(WNG),WNG,nw);
%             qd=QBild(W1D,W6D,W7D,Pt(vd)',Uf(WonD(:,3),ft+1),WonD(:,1),dt,pw(WonD(:,3)),WonD(:,3),nw);
%             q=qm+qc+qg+qd;
%             flag_pwq=0;
%            % [flag_pwq,Pw(:,ft+1),Qz(:,ft+1),Qf]=Chek_bond2(pw,Pt(Won),Uf(WonM,ft+1),Qf,PwQC_bnd,q/dt);
%             Pt0=Pt;
%         end
%         Pi(:,1)=Pt(1:na+nc+ng+nd+nb);    
%     end
%     
%     Pi0=Pi;
%         
%     Pw(Qf~=0,ft+1)=Pt(na+nc+ng+nd+nb+1:end);
%     Pwt(:,t+1)=Pw(:,ft+1);
%     %% Водонасыщенность
%     
%     %[SCw,SCp,NDT(t)]=Sat_fast(SCw,SCp,RC,TC,TG,TA2C,TA2G,TM,Pi(:,1),PR,ndt,Won,Wf,...
%     %    Uf(:,t),dt,dVCG,Pw(:,t),WonG,CpW(:,t),WonC,Nl,b2gm,GYData.GY_Pz);
%     
% % %     qm=QBild(W1,W6,W7,Pi(1:na,1),Uf(WonM,ft+1),Won,dt,Pw(WonM,ft+1),WonM,nw);
% % %     qc=QBild(W1C,W6C,W7C,Pi(vc,1),Uf(WonC(:,3),ft+1),WonC(:,1),dt,Pw(WonC(:,3),ft+1),WonC(:,3),nw);
% % %     qg=QBild(W1G,W6G,W7G,Pi(vg,1),Uf(WNG,ft+1),WonG(:,1),dt,Pw(WNG,ft+1),WNG,nw);
% % %     qd=QBild(W1D,W6D,W7D,Pi(vd,1),Uf(WonD(:,3),ft+1),WonD(:,1),dt,Pw(WonD(:,3),ft+1),WonD(:,3),nw);
% % %     q=qm+qc+qg+qd;
%     Qz1=q(:,1)+q(:,2);
%     % Qzm1=qm(:,1)+qm(:,2);
%     
%     Qf=Qz1;
%     
%   if (isempty(RC.Cr)==0 || isempty(RC.Gr)==0) && ft==0
%     [ndt,~,~]=vibor_t(ndt,0,Pi(vc),Pi(vg),Pw(:,ft+1),PR,RC,dt,Sw0([vc,vg]),Sw([vc,vg]),...
%         CL(RC.Cr2+(RC.Cc2-1)*nc),GL(RC.Gr2+(RC.Gc2-1)*ng),dVCG(vc),dVCG(vg),W1C,W1G,WonC,WonG,ft+1,0);
%   end;
%     Sw0=Sw;
%     [Sw,Cp,NDT(t),Q1,Q2,Qm1,Qd1,dSS(t),ndt]=Sat_fast_2(Sw,Cp,RC,TC,TG,TA2C,TA2G,TA2D,TD2C,TD2G,Pi(:,1),PR,ndt,Won,...
%         Uf(:,ft+1),dt,dVCG,Pw(:,ft+1),WonG,CpW(:,ft+1),WonC,Nl,CR_rc,Qz1,Qf,Pi0,TL,W1,TW,W6,TP,...
%         W7,L,DATA.Lc,Lg,Ke,Cws,Cwp,BLGY_GIM,Qz(:,ft+1),WonM,nw,b1gm,b1gd,GYData,Clp,ka1,...
%         W1D,W6D,W7D,A2BW,A2BP,D2BW,D2BP,DL,DW,DP,WonD,A2DW,A2DP,A2DL,BB,A2BL,D2BL,b1gb,Grw,dZ,Mp,Bwo,P0);
%     %
%     %  [SCw,SCp,NDT(:,t),Q1,Q2,Qm1,dSS(t)]=Sat_fast_1(SCw,SCp,RC,TC,TG,TA2C,TA2G,Pi(:,1),PR,ndt,Won,...
%     %      Uf(:,ft+1),dt,dVCG,Pw(:,ft+1),WonG,CpW(:,ft+1),WonC,Nl,CR_ind,Qz1,Qf,Pi0,TL,W1,TW,W6,TP,...
%     %      W7,L,DATA.Lc,DATA.Lg,Ke,Cws,Cwp,BLGY_GIM,Qz(:,ft+1),WonM,nw,b1gm,Clp);
%     
%     MSw(:,1)=Sw(va);
%     CSw(:,t+1)=Sw(vc);
%     GSw(:,t+1)=Sw(vg);
%     DSw(:,1)=Sw(vd);
%     
%     MCp(:,1)=Cp(va);
%     CCp(:,t+1)=Cp(vc);
%     GCp(:,t+1)=Cp(vg);
%     DCp(:,1)=Cp(vd);
%     
%     %% Дебиты
%     Qm(:,:,t+1)=QBild(W1,W6.*A(Won(:,1)),W7,Pi(va,1),Uf(WonM,ft+1),Won(:,1),dt,Pw(WonM,ft+1),WonM,nw);
%     %Qm(CR_rc.wn,:,t+1)=Qm1(CR_rc.wn,:);
%     Qm(CR_rc(1,1).won(:,3),:,t+1)=Qm1(CR_rc(1,1).won(:,3),:);
%     Qc(:,:,t+1)=Q1;%QBild(W1C,W6C,W7C,Pi(na+1:na+nc,1),Uf(WonC(:,3),t+1),WonC(:,1),dt,Pw(WonC(:,3),t+1));
%     Qg(:,:,t+1)=Q2;%QBild(W1G,W6G,W7G,Pi(na+nc+1:end,1),Uf(WNG,t+1),WonG(:,1),dt,Pw(WNG,t+1));
%     Qd(:,:,t+1)=QBild(W1D,W6D.*A(WonD(:,1)),W7D,Pi(vd,1),Uf(WonD(:,3),ft+1),WonD(:,1),dt,Pw(WonD(:,3),ft+1),WonD(:,3),nw);
%     Qd(CR_rc(1,2).won(:,3),:,t+1)=Qd1(CR_rc(1,2).won(:,3),:);
%     
%     PpW(WonM,t+1)=Pi(Won(:,1));
%     %% пїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ
%     %    Ti(:,t+1)=Termal(Ti(:,t),Pi(:,t+1),SCw,Mp,Mc,Mg,NTG,LM,LC,LG,LA2C,LA2G,PR,RC,dVCG,dt,Won,WonM,WonV,WonG,...
%     %       TeW(:,t+1),Qm(:,:,t+1),Qc(:,:,t+1),Qg(:,:,t+1),C2GL,TL,CL,GL,TW,CW,GW,A2CL,A2CW,A2GL,A2GW,Wf,NDT(t),t,T,S,BZ);
%     
%     Ti(:,1)=Ti(:,1);
%     
%     if mod(t,zapt)==0
%         if Ta/dt<=100
%             zapt=1;
%         elseif Ta/dt<=1000
%             zapt=10;
%         elseif Ta/dt<=1e4
%             zapt=100;
%         else
%             zapt=1000;
%         end;
%         
%         j=j+1;
%         Pj(:,j)=Pi;
%         Swj(:,j)=MSw;
%         DSwj(:,j)=DSw;
%         CSwj(:,j)=Sw(vc);
%         MCpj(:,j)=MCp;
%         Tj(:,j)=Ti;
%         waitbar(st/Ta)
%     end;
%     sQo=sum(Qm(:,2,t+1)-Qm(:,3,t+1)+Qm(:,1,t+1))+sum(Qc(:,2,t+1)-Qc(:,3,t+1)+Qc(:,1,t+1))...
%         +sum(Qg(:,2,t+1)-Qg(:,3,t+1)+Qg(:,1,t+1))+sum(Qd(:,2,t+1)-Qd(:,3,t+1)+Qd(:,1,t+1));
%     %sum(sQo(:))
%     % dQ(t)=sum(Sw0.*[dV;dVC;dVG])-sum([Sw;Cw(:,t+1);Gw(:,t+1)].*[dV;dVC;dVG])-sum(sQo(:));
%     %%prob.progress;
%     c_lik=1-Qm(:,3,t+1)./Qm(:,2,t+1);
%     c_lik(isnan(c_lik)==1)=0;
%     Uf(:,ft+1:end)=Uf(:,ft+1:end).*repmat(c_lik<PwQC_bnd(:,7),1,size(Uf(:,ft+1:end),2));
%     qin=-Qm(:,1,t+1)/dt;
%     qo=Qm(:,3,t+1)/dt;
%     qo(Uf(:,ft+1)==-1)=inf;
%     Qmin=repmat(PwQC_bnd(:,8),1,1);
%     Uf(:,ft+1:end)=Uf(:,ft+1:end).*repmat(qo>=Qmin,1,size(Uf(:,ft+1:end),2));
%     
%     st=st+dt;
%     dt1(t+1)=dt;
%     dt=vibor_t2(dtt,Pi,RC,dVCG,TL,W1,Won(:,1),Pw(:,ft+1),na,PR,st,Ta,Sw,Sw0,dt,Nl,WonM,va,vd,DL,W1D,WonD,nd,dV1,dV2);
%     t_flag=st~=Ta;
%     
% end;
% 
% j=j+1;
% Pj(:,j)=Pi;
% Swj(:,j)=MSw;
% DSwj(:,j)=DSw;
% CSwj(:,j)=Sw(vc);
% MCpj(:,j)=MCp;
% Tj(:,j)=Ti;
% waitbar(st/Ta)
% [Q,Pw,PpW,dtz]=Q2Sut(Qm,Qc,Qg,Qd,Pwt(:,1:t+1),PpW(:,1:t+1),dt1,Ta);
% % GY_Pxy(p)=GY_Pxy;
% % save('GY_Pxy.mat','GY_Pxy')
% %  Bnd_xy(p)=DATA.BndXY;
% %  GY_Pxy(Bnd_xy==1)'
% toc