function [Pj,Swj,Tj,MCpj,p,Q,Pw,PpW,CSw,NDT,Uf,dt1,dV0,ka,dtz,DSwj]=SimT_MKT(PR,C,A2C,GData,BB,A2B,DData,dVC,dVB,DATA,WData,GYData,fll,CR_GRUP)
tic

KX=DATA.gKX;
KY=DATA.gKY;
KZ=DATA.gKZ;
Mp=DATA.gMp;
P=DATA.gP;
MSw=DATA.gSw;
MCp=DATA.gCp;
XY=DATA.XY;
BND=DATA.BND;
BndXY=DATA.BndXY;
BndZ=DATA.BndZ;
H=DATA.gH;
Z=DATA.gZ;
ka=DATA.ka;

[G,A2G,dVG,WonG,Lg,Mp_g]=Rasp(GData,3);   %Распаковка структуры для гор. трещ.
[D,A2D,dVD,WonD,Ld,Mp_d]=Rasp(DData,4);   %Распаковка структуры для двойной среды

Won=DATA.Won;
WNG=WonG(:,3);
WonC=DATA.WonV;

Pw=WData.Pw;
Uf=WData.Uf;
CpW=WData.CpW;
TeW=WData.TeW;
Qz=WData.Qz;
PwQC_bnd=WData.PwQC_bnd;

Nl=PR.Nl;  as=PR.as; aw=PR.aw; ts=PR.ts; tw=PR.tw; mu=PR.mu; mup=PR.mup;
Ta=PR.Ta;  dtt=PR.dt; ndt=PR.ndt; zc=PR.zc; Bo=PR.Bo;  kms=PR.kms;  Ro=PR.Ro;
nw=size(Qz,1);

Qm=zeros(nw,5,size(Uf,2));
Qc=zeros(nw,5,size(Uf,2));
Qg=zeros(nw,5,size(Uf,2));
Qd=zeros(nw,5,size(Uf,2));

qm=zeros(size(Uf,1),5);
qc=zeros(size(Uf,1),5);
qg=zeros(size(Uf,1),5);
qd=zeros(size(Uf,1),5);

%KX(:)=mean(KX(:));

[A]=MR_Prop_Bond(XY,Nl,BND);
[L,B,S,H1,HV]=Geome3_1(A,XY,Z,H);
%ka(sum(A)==-1)=0;
[r,c]=find(A(1:size(XY,1),1:size(XY,1))==1);

Wf=KWell(KX,H,S,L,B,Won,r,c,WData.Doly,WData.SDoly,WData.r0,XY,nw,Nl);
XY=repmat(XY,Nl,1);
ka1=ka(Won);
Won1=Won;
[A,L,B,S,H1,HV,XY,KX,KY,KZ,Mp,MSw,H,Z,P,MCp,Won,Wf]=YdalActcell(A,L,B,S,H1,HV,XY,KX,KY,KZ,Mp,MSw,H,Z,P,MCp,Won,Wf,ka);
[A,L,S,B,H1,HV,K,XY,Mp,MSw,H,Z,P,MCp,T,NTG,p,rz,cz,BXY,BZ,dH,NL,NamXY,GYData]=PereYpor(A,L,S,B,H1,HV,KX,KY,KZ,Mp,MSw,...
    XY,H,Z,P,MCp,DATA,GYData,ka);

[r,c]=find(L);
nw1=size(p,2);
nwc=size(A2C,1);
raz=nw1-nwc;

A2C(raz+1:end+raz,:)=A2C;
A2C(1:raz,:)=0;
A2C=A2C(p,:);
A2G=A2G(p,:);
A2D=A2D(p,:);
A2B=A2B(p,:);

for i=1:size(Won,1)
    Won(i)=find(Won(i)==p);
end;

na=size(XY,1);  nc=size(C,1);  ng=size(G,1);  nd=size(D,1); nb=size(BB,1); nw=size(Qz,1);
va=1:na;
vc=na+1:na+nc;
vg=na+nc+1:na+nc+ng;
vd=na+nc+ng+1:na+nc+ng+nd;
vb=na+nc+ng+nd+1:na+nc+ng+nd+nb;

C2G=sparse(nc,ng);     C2GL=C2G;
C2B=sparse(nc,nb);     C2BL=C2B;
G2B=sparse(ng,nb);     G2BL=G2B;
D2B=sparse(nd,nb);
b1wb=sparse(nb,1);

C2D=DData.C2D;     D2CL=C2D;
G2D=DData.G2D;     D2GL=G2D;

WonM=repmat(1:nw,1,Nl)';
WW=WonM~=0;
WonM=WonM(ka1==1);
WW(ka1==0)=0;
WW=reshape(WW,Nl,nw);
WW=sum(WW,1)~=0;

%figure(98),subplot(2,4,7),spy([A,A2C,A2G;A2C',C,C2G;A2G',C2G',G]);

[Ke,Ke_gy,dV]=KH2Mat(K,HV,Mp,S,r,c,rz,cz,GYData.GY_Kz,GYData.GY_Kxy,BndXY(p),DData);

Mp_c=ones(nc,1);
Mp_b=100*ones(nb,1);

CCp(:,1)=zeros(nc,1);
GCp(:,1)=zeros(ng,1);
DCp(:,1)=zeros(nd,1);
BCp(:,1)=zeros(nb,1);

dVCG=[dV;dVC;dVG;dVD.*Mp_d;dVB];
dV0=[dV;dVC;dVG;dVD.*Mp_d];

[TM,TC,TG,TD,TA2C,TA2G,TA2D,TD2C,TD2G,RC,Txyz_GY_A,Txyz_GY_D,dV1,dV2,TTM]=Pre_fast(A,C,G,D,A2C,A2G,A2D,C2G,C2D,G2D,...
    Ke,L,B,S,H1,K(:,1),DData.Kd,Ke_gy,BndXY(p),BndZ(p),nb,dVCG);
[dZ]=dZ_bild(Z,RC,PR.g,PR.Ro);

vad=RC.ADr;
CSw(:,1)=MSw(RC.ACr,1);
GSw(:,1)=MSw(RC.AGr,1);
DSw(:,1)=MSw(RC.ADr,1);
BSw(:,1)=ones(nb,1);

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

j=0;
Pwt=zeros(size(Pw));
PpW=zeros(size(Uf));

% [CR_rc]=Pre_Crack1(RC,na,TM,A2C,A2G,Wf,Won,WonM,nw);
[CR_rc]=Pre_Crack(RC,na,nd,TTM,TD,A2C,A2G,C2D,G2D,A2D,Wf,Won,WonM,WonD,dZ);
%[CR_ind]=Pre_Crack_p(RC,na,TM,Wf,Won,WonM,nw,CR_GRUP,C,G,A2C,A2G,K(:,1),WonC,WonG);
[CR]=SS_ind(RC,na);

st=0;
t=0;
t_flag=1;
dt=dtt;

Bwo(:,1:2)=Bo(1)*ones(na+nc+ng+nd+nb,2);
Bwo(:,3:4)=Bo(2)*ones(na+nc+ng+nd+nb,2);

Mp=repmat([Mp.*ones(na,1);Mp_c;Mp_g;Mp_d;Mp_b],1,2);

Sw=[MSw(:,1);CSw(:,1);GSw(:,1);DSw(:,1);BSw(:,1)];
Cp=[MCp(:,1);CCp(:,1);GCp(:,1);DCp(:,1);BCp(:,1)];

V0=sum(dVCG.*(1-Sw));
P0=Pi;
Pt0=[P0;Pw(:,1)]';
dtt1=dtt;
Wf0=Wf;
GY_Pxy=GYData.GY_Pxy;
GY_Pz=GYData.GY_Pz;
GY_Pxy2=GYData.P0;
zapt=1;

while t_flag==1
    %for t=1:Ta-1
    t=t+1;
    
    if t==1 && dtt==0
        dt=1;
        dtt=1;
    else
        dtt=dtt1;
    end;
    
    ft=floor(st);
    % if ~isempty(find(CpW(:,t)~=0)) fp=1; else fp=0; end;
    fp=1;
    
    Qf=Qz(:,ft+1);
    Qf(WW==0)=0;
    
    %     [C,G,A2C,A2G,WonC(:,2),WonG(:,2)]=temp_CG(fll,C2,G2,A2C2,A2G2,WonC1,WonG1,t);
    %     [TC,TG,TA2C,TA2G,WonC(:,2),WonG(:,2)]=temp_CG(fll,TC1,TG1,TA2C1,TA2G1,WonC1,WonG1,t);
    
    %     [C,TC]=GeoMex(C2,TC1,Pi(na+1:na+nc,1),1);
    %     [G,TG]=GeoMex(G2,TG1,Pi(na+nc+1:na+nc+ng,1),1);
    
    kfw(1:na,1)=Sat_cal(MSw,1,1,as,aw); %water
    kfo(1:na,1)=Sat_cal(MSw,2,1,as,aw); %oil
    
    kfw(na+1:na+nc+ng,1)=Sat_cal([CSw(:,t);GSw(:,t)],1,1,ts,tw); %water
    kfo(na+1:na+nc+ng,1)=Sat_cal([CSw(:,t);GSw(:,t)],2,1,ts,tw); %oil
    
    [b1gm,b1gc,b1gg,b1gd,b1gb]=GY_bild(GYData,Pi([va,vd],1),Sw([va,vd],1),Cp([va,vd],1),RC,Txyz_GY_A,Txyz_GY_D,PR);
    
    kj=0;
    flag_gim=1;
    flag_pwq=1;
    
    while flag_gim==1 && kj<10
        kj=kj+1;
        
        [Clp,Cwp,Cws,A,Bwo,Mp]=SGim(dVCG./Mp(:,1),Sw,Mp,zc,Bwo,Pi,1,P0,va,vc,vg,vd,vb,dt);
        
        [TL,TW,TP]=Potok_MKT(TTM,Pi(va,1),kfw(va),kfo(va),MCp(:,1),mu,RC.Arc2,mup,fp,kms(1),L,Ke,Ro,A(va),dZ(1,:));
        [CL,CW,~]=Potok_Tube(TC,Pi(vc,1),CSw(:,t),CCp(:,t),PR,mup,fp,kms(2),DATA.Lc,RC.Cr2,RC.Cc2,nc,A(vc),dZ(2,:));
        [GL,GW,~]=Potok_Tube(TG,Pi(vg,1),GSw(:,t),GCp(:,t),PR,mup,fp,kms(3),Lg,RC.Gr2,RC.Gc2,ng,A(vg),dZ(3,:));
        [DL,DW,DP]=Potok_Tube(TD,Pi(vd,1),DSw(:,1),DCp(:,1),PR,mup,fp,kms(4),Ld,RC.Dr2,RC.Dc2,nd,A(vd),dZ(4,:));
        
        [Gr,Grw]=Gravity(TL,TW,CL,CW,GL,GW,DL,DW,[],A,dZ);
        
        [A2CL,~,~]=Obmen_T2M(A2C,Pi(va,1),Pi(vc,1),MSw(:,1),CSw(:,t),K(:,1),PR,MCp(:,1),CCp(:,t));
        [A2GL,~,~]=Obmen_T2M(A2G,Pi(va,1),Pi(vg,1),MSw(:,1),GSw(:,t),K(:,1),PR,MCp(:,1),GCp(:,t));
        
        [A2DL,A2DW,A2DP]=Obmen_T2M(A2D,Pi(va,1),Pi(vd,1),MSw(:,1),DSw(:,1),K(:,1),PR,MCp(:,1),DCp(:,1));
        [A2BL,A2BW,A2BP]=Obmen_T2M(A2B,Pi(va,1),Pi(vb,1),MSw(:,1),BSw(:,1),ones(na,1),PR,MCp(:,1),BCp(:,1));
        [D2BL,D2BW,D2BP]=Obmen_T2M(D2B,Pi(vd,1),Pi(vb,1),DSw(:,1),BSw(:,1),K(:,1),PR,DCp(:,1),BCp(:,1));
        
        [D2CL,~,~]=Obmen_T2M(C2D,Pi(vd,1),Pi(vc,1),DSw(:,1),CSw(:,t),K(:,1),PR,DCp(:,1),CCp(:,t));
        [D2GL,~,~]=Obmen_T2M(G2D,Pi(vd,1),Pi(vg,1),DSw(:,1),GSw(:,t),K(:,1),PR,DCp(:,1),GCp(:,t));
        
        %Wf=Wf0.*(1-0.0001*(P0(Won)-Pi(Won))).^3;
        [W1,W6,W7]=Well_MKT(Wf,Won,Uf(WonM,ft+1),MSw(:,1),MCp(:,1),aw,as,mu,mup,CpW(WonM,ft+1),A(va));
        [W1C,W6C,W7C]=Well_MKT(WonC(:,2),WonC(:,1),Uf(WonC(:,3),ft+1),CSw(:,t),CCp(:,t),tw,ts,mu,mup,CpW(WonC(:,3),ft+1),A(vc));
        [W1G,W6G,W7G]=Well_MKT(WonG(:,2),WonG(:,1),Uf(WNG,ft+1),GSw(:,t),GCp(:,t),tw,ts,mu,mup,CpW(WNG,ft+1),A(vg));
        [W1D,W6D,W7D]=Well_MKT(WonD(:,2),WonD(:,1),Uf(WonD(:,3),ft+1),DSw(:,1),DCp(:,1),tw,ts,mu,mup,CpW(WonD(:,3),ft+1),A(vd));
        
        A1=TL-sparse(Won,Won,W1,na,na)-sparse(1:na,1:na,sum(A2CL,2)+sum(A2GL,2)+sum(A2BL,2)+sum(A2DL,2)+Clp(va)+sum(b1gm(:,1:2),2),na,na);  %Матрица коэф. для пор
        C1=CL-sparse(1:nc,1:nc,sum(A2CL,1)+sum(D2CL,1)+Clp(vc)',nc,nc)-sparse(WonC(:,1),WonC(:,1),W1C,nc,nc);                       %Матрица коэф. для вертикальных трещ.
        G1=GL-sparse(1:ng,1:ng,sum(A2GL,1)+sum(D2GL,1)+Clp(vg)',ng,ng)-sparse(WonG(:,1),WonG(:,1),W1G,ng,ng);                       %Матрица коэф. для гориз. трещ.
        D1=DL-sparse(WonD(:,1),WonD(:,1),W1D,nd,nd)-sparse(1:nd,1:nd,sum(A2DL,1)'+sum(D2BL,2)+sum(D2CL,2)+sum(D2GL,2)+Clp(vd)+sum(b1gd(:,1:2),2),nd,nd);                       %Матрица коэф. для двойной пор.
        B1=BB-sparse(1:nb,1:nb,sum(A2BL,1)+sum(D2BL,1)+Clp(vb)'+b1gb',nb,nb);                                                             %Матрица коэф. для границ
        
        W2M=sparse(WonM,Won,W1,nw,na);
        W2C=sparse(WonC(:,3),WonC(:,1),W1C,nw,nc);
        W2G=sparse(WNG,WonG(:,1),W1G,nw,ng);
        W2D=sparse(WonD(:,3),WonD(:,1),W1D,nw,nd);
        W2B=sparse(nw,nb);

        WM1=[W2M,W2C,W2G,W2D,W2B];
        WM2=WM1';
        W3vec=sparse(WonM,1,W1,nw,1)+sparse(WonC(:,3),1,W1C,nw,1)+sparse(WNG,1,W1G,nw,1)+sparse(WonD(:,3),1,W1D,nw,1);
        WM3=-sparse(1:nw,1:nw,W3vec,nw,nw);
        
        AM=[A1,   A2CL, A2GL, A2DL, A2BL;
            A2CL',  C1,  C2GL, D2CL', C2BL;
            A2GL', C2GL', G1,  D2GL', G2BL;
            A2DL', D2CL, D2GL,  D1,  D2BL;
            A2BL', C2BL',G2BL',D2BL', B1];
            
        while  flag_pwq==1
           
            PwNl=repmat(Pw(:,ft+1),Nl,1);
            % PwNl=PwNl(ka1==1);
            
            b1wm=sparse(Won,ones(1,size(Won,1)),-W1.*PwNl(WonM),na,1);
            b1wc=sparse(WonC(:,1),ones(1,size(WonC,1)),-W1C.*Pw(WonC(:,3),ft+1),nc,1);
            b1wg=sparse(WonG(:,1),ones(1,size(WonG,1)),-W1G.*Pw(WNG,ft+1),ng,1);
            b1wd=sparse(WonD(:,1),ones(1,size(WonD,1)),-W1D.*PwNl(WonD(:,3)),nd,1);
            
            b1wm=b1wm.*(sum(W2M(Qf==0,:),1)~=0)';
            b1wc=b1wc.*(sum(W2C(Qf==0,:),1)~=0)';
            b1wg=b1wg.*(sum(W2G(Qf==0,:),1)~=0)';
            b1wd=b1wd.*(sum(W2D(Qf==0,:),1)~=0)';
            %  b1wm

            BM=[b1wm;b1wc;b1wg;b1wd;b1wb]+[-b1gm(:,2).*GY_Pz-b1gm(:,1).*GY_Pxy;b1gc;b1gg;...
                [-b1gd(:,2).*GY_Pz(vad)-b1gd(:,1).*GY_Pxy(vad)];-b1gb.*GY_Pxy2]-(Clp.*Pi)-Gr;
            
            BLGY_GIM=[[-b1gm(:,2).*GY_Pz-b1gm(:,1).*GY_Pxy];[-b1gd(:,2).*GY_Pz(vad)-b1gd(:,1).*GY_Pxy(vad)];...
                -b1gb.*GY_Pxy2]-(Clp([va,vd,vb]).*Pi([va,vd,vb]))-Gr([va,vd,vb]);
            
            W2M1=WM1(Qf~=0,:);
            W2M2=WM2(:,Qf~=0);
            W2M3=WM3(Qf~=0,Qf~=0);
            
            Pt=[BM',Qz(Qf~=0,ft+1)']/[AM,W2M2;W2M1,W2M3];
            
            flag_gim=sum(abs(Pt(1:na+nc+ng+nd+nb)-Pt0(1:na+nc+ng+nd+nb))./Pt(1:na+nc+ng+nd+nb)>=1e-6)~=0;
            flag_gim=flag_gim*(sum(zc==0)==0);
            
            pw=Pw(:,ft+1);
            pw(Qf~=0)=Pt(na+nc+ng+nd+nb+1:end);
            
            %[flag_pwq,Pw(:,ft+1),Qz(:,ft+1),Qf]=Chek_bond(pw,Pt(Won),W1,Uf(WonM,ft+1),Qf,PwQC_bnd);
            
            qm=QBild(W1,W6,W7,Pt(va)',Uf(WonM,ft+1),Won,dt,pw(WonM),WonM,nw);
            qc=QBild(W1C,W6C,W7C,Pt(vc)',Uf(WonC(:,3),ft+1),WonC(:,1),dt,pw(WonC(:,3)),WonC(:,3),nw);
            qg=QBild(W1G,W6G,W7G,Pt(vg)',Uf(WNG,ft+1),WonG(:,1),dt,pw(WNG),WNG,nw);
            qd=QBild(W1D,W6D,W7D,Pt(vd)',Uf(WonD(:,3),ft+1),WonD(:,1),dt,pw(WonD(:,3)),WonD(:,3),nw);
            q=qm+qc+qg+qd;
            
            [flag_pwq,Pw(:,ft+1),Qz(:,ft+1),Qf]=Chek_bond2(pw,Pt(Won),Uf(WonM,ft+1),Qf,PwQC_bnd,q/dt);
            Pt0=Pt;
        end
    end
    
    Pi0=Pi;
    Pi(:,1)=Pt(1:na+nc+ng+nd+nb);
    
    Pw(Qf~=0,ft+1)=Pt(na+nc+ng+nd+nb+1:end);
    Pwt(:,t+1)=Pw(:,ft+1);
    %% Водонасыщенность
    
    %[SCw,SCp,NDT(t)]=Sat_fast(SCw,SCp,RC,TC,TG,TA2C,TA2G,TM,Pi(:,1),PR,ndt,Won,Wf,...
    %    Uf(:,t),dt,dVCG,Pw(:,t),WonG,CpW(:,t),WonC,Nl,b2gm,GYData.GY_Pz);
    
% %     qm=QBild(W1,W6,W7,Pi(1:na,1),Uf(WonM,ft+1),Won,dt,Pw(WonM,ft+1),WonM,nw);
% %     qc=QBild(W1C,W6C,W7C,Pi(vc,1),Uf(WonC(:,3),ft+1),WonC(:,1),dt,Pw(WonC(:,3),ft+1),WonC(:,3),nw);
% %     qg=QBild(W1G,W6G,W7G,Pi(vg,1),Uf(WNG,ft+1),WonG(:,1),dt,Pw(WNG,ft+1),WNG,nw);
% %     qd=QBild(W1D,W6D,W7D,Pi(vd,1),Uf(WonD(:,3),ft+1),WonD(:,1),dt,Pw(WonD(:,3),ft+1),WonD(:,3),nw);
% %     q=qm+qc+qg+qd;
    Qz1=q(:,1)+q(:,2);
    
    Qzm1=qm(:,1)+qm(:,2);
    
    Qf=Qz1;
    Sw0=Sw;
    
    [Sw,Cp,NDT(t),Q1,Q2,Qm1,Qd1,dSS(t)]=Sat_fast_2(Sw,Cp,RC,TC,TG,TA2C,TA2G,TA2D,TD2C,TD2G,Pi(:,1),PR,ndt,Won,...
        Uf(:,ft+1),dt,dVCG,Pw(:,ft+1),WonG,CpW(:,ft+1),WonC,Nl,CR_rc,Qz1,Qf,Pi0,TL,W1,TW,W6,TP,...
        W7,L,DATA.Lc,Lg,Ke,Cws,Cwp,BLGY_GIM,Qz(:,ft+1),WonM,nw,b1gm,b1gd,GYData,Clp,ka1,...
        W1D,W6D,W7D,A2BW,A2BP,D2BW,D2BP,DL,DW,DP,WonD,A2DW,A2DP,A2DL,BB,A2BL,D2BL,b1gb,Grw,dZ,Mp,Bwo,P0);
    %
    %  [SCw,SCp,NDT(:,t),Q1,Q2,Qm1,dSS(t)]=Sat_fast_1(SCw,SCp,RC,TC,TG,TA2C,TA2G,Pi(:,1),PR,ndt,Won,...
    %      Uf(:,ft+1),dt,dVCG,Pw(:,ft+1),WonG,CpW(:,ft+1),WonC,Nl,CR_ind,Qz1,Qf,Pi0,TL,W1,TW,W6,TP,...
    %      W7,L,DATA.Lc,DATA.Lg,Ke,Cws,Cwp,BLGY_GIM,Qz(:,ft+1),WonM,nw,b1gm,Clp);
    
    MSw(:,1)=Sw(va);
    CSw(:,t+1)=Sw(vc);
    GSw(:,t+1)=Sw(vg);
    DSw(:,1)=Sw(vd);
    
    MCp(:,1)=Cp(va);
    CCp(:,t+1)=Cp(vc);
    GCp(:,t+1)=Cp(vg);
    DCp(:,1)=Cp(vd);
    
    %% Дебиты
    
    Qm(:,:,t+1)=QBild(W1,W6.*A(Won),W7,Pi(va,1),Uf(WonM,ft+1),Won,dt,Pw(WonM,ft+1),WonM,nw);
    %Qm(CR_rc.wn,:,t+1)=Qm1(CR_rc.wn,:);
    Qm(CR_rc(1,1).won(:,3),:,t+1)=Qm1(CR_rc(1,1).won(:,3),:);
    Qc(:,:,t+1)=Q1;%QBild(W1C,W6C,W7C,Pi(na+1:na+nc,1),Uf(WonC(:,3),t+1),WonC(:,1),dt,Pw(WonC(:,3),t+1));
    Qg(:,:,t+1)=Q2;%QBild(W1G,W6G,W7G,Pi(na+nc+1:end,1),Uf(WNG,t+1),WonG(:,1),dt,Pw(WNG,t+1));
    Qd(:,:,t+1)=QBild(W1D,W6D.*A(WonD(:,1)),W7D,Pi(vd,1),Uf(WonD(:,3),ft+1),WonD(:,1),dt,Pw(WonD(:,3),ft+1),WonD(:,3),nw);
    Qd(CR_rc(1,2).won(:,3),:,t+1)=Qd1(CR_rc(1,2).won(:,3),:);
    
    PpW(WonM,t+1)=Pi(Won);
    %% пїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ
    %    Ti(:,t+1)=Termal(Ti(:,t),Pi(:,t+1),SCw,Mp,Mc,Mg,NTG,LM,LC,LG,LA2C,LA2G,PR,RC,dVCG,dt,Won,WonM,WonV,WonG,...
    %       TeW(:,t+1),Qm(:,:,t+1),Qc(:,:,t+1),Qg(:,:,t+1),C2GL,TL,CL,GL,TW,CW,GW,A2CL,A2CW,A2GL,A2GW,Wf,NDT(t),t,T,S,BZ);
    
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
    
    dt=vibor_t2(dtt,Pi,RC,dVCG,TL,W1,Won,Pw(:,ft+1),na,PR,st,Ta,Sw,Sw0,dt,Nl,WonM,va,vd,DL,W1D,WonD,nd,dV1,dV2);
    st=st+dt;
    t_flag=st~=Ta;
    dt1(t)=dt;
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