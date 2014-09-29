function [XY,KX,Z,Pj,Swj,Tj,MCpj,p,Q,Pw,PpW,Cw,NDT,dQ,dSS,dt1,V0]=SimT_MKT(PR,C,A2C,G,A2G,dVC,dVG,DATA,WData,GYData,fll)
tic

KX=DATA.gKX;
KY=DATA.gKY;
KZ=DATA.gKZ;
Mp=DATA.gMp;
P=DATA.gP;
Sw=DATA.gSw;
MCp=DATA.gCp;
XY=DATA.XY;
H=DATA.gH;
Z=DATA.gZ;


Won=DATA.Won;
WonG=DATA.WonG;
WNG=DATA.WonG(:,3);
WonC=DATA.WonV;

Pw=WData.Pw;
Uf=WData.Uf;
CpW=WData.CpW;
TeW=WData.TeW;
Qz=WData.Qz;

Nl=PR.Nl;  as=PR.as; aw=PR.aw; ts=PR.ts; tw=PR.tw; mu=PR.mu; mup=PR.mup;
Ta=PR.Ta;  dtt=PR.dt; ndt=PR.ndt; zc=PR.zc; Bo=PR.Bo;  kms=PR.kms;  Ro=PR.Ro;

Qm=zeros(size(Uf,1),5,size(Uf,2));
Qc=zeros(size(WonC,1),5,size(Uf,2));
Qg=zeros(size(WNG,1),5,size(Uf,2));

qm=zeros(size(Uf,1),5);
qc=zeros(size(Uf,1),5);
qg=zeros(size(Uf,1),5);

%KX(:)=mean(KX(:));

[A]=MR_Prop(XY,Nl);

[L,B,S,H1]=Geome3_1(A,XY,Z,H);

[r,c]=find(A==1);

Wf=KWell(KX,H,S,L,B,Won,r,c,WData.Doly,WData.r0,XY);
XY=repmat(XY,Nl,1);

[A,L,S,B,H1,K,XY,Mp,Sw,H,Z,P,MCp,T,NTG,p,rz,cz,BXY,BZ,dH,NL,NamXY,GYData]=PereYpor(A,L,S,B,H1,KX,KY,KZ,Mp,Sw,XY,H,Z,P,MCp,DATA,GYData);
[r,c]=find(L);

nw=size(p,2);
nwc=size(A2C,1);
raz=nw-nwc;
A2C(raz+1:end+raz,:)=A2C;
A2C(1:raz,:)=0;
A2C=A2C(p,:);
A2G=A2G(p,:);

for i=1:size(Won,1)
  Won(i)=find(Won(i)==p);
end;

na=size(XY,1);  nc=size(C,1);  ng=size(G,1);    nw=size(Qz,1);  nd=0;
C2G=sparse(nc,ng);
C2GL=C2G;
WonM=repmat(1:nw,1,Nl);

figure(98),subplot(2,4,7),spy([A,A2C,A2G;A2C',C,C2G;A2G',C2G',G]);

[Ke,Ke_gy,dV]=KH2Mat(K,H,Mp,S,r,c,rz,cz,GYData.GY_Kz,BZ); 

Mc=ones(nc,1);
Mg=ones(ng,1);

CCp(:,1)=zeros(nc,1);
GCp(:,1)=zeros(ng,1);

dVD=[];
dVCG=[dV;dVC;dVG;dVD];

[TM,TC,TG,TA2C,TA2G,RC,T_GY]=Pre_fast(A,C,A2C,A2G,C2G,G,Ke,L,B,H1,K(:,1),Ke_gy,BZ,rz,cz);
Cw(:,1)=Sw(RC.ACr,1);
Gw(:,1)=Sw(RC.AGr,1);
Pi(1:na,1)=P;
Pi(na+1:na+nc,1)=P(RC.ACr);
Pi(na+nc+1:na+nc+ng,1)=P(RC.AGr);

Ti(1:na,1)=T;
Ti(na+1:na+nc,1)=T(RC.ACr);
Ti(na+nc+1:na+nc+ng,1)=T(RC.AGr);

j=0;
PpW=zeros(size(Uf));

[CR_rc]=Pre_Crack(RC,na,TM,A2C,A2G,Wf,Won,WonM,nw);
[CR]=SS_ind(RC,na);


C2=C;
G2=G;
A2C2=A2C;
A2G2=A2G;

TC1=TC;
TG1=TG;
TA2C1=TA2C;
TA2G1=TA2G;
%WonC(:,2)=0;
WonC1=WonC(:,2);
WonG1=WonG(:,2);

st=0;
t=0;
t_flag=1;
Sw2=Sw;
dt=dtt;
va=1:na;
vc=na+1:na+nc;
vg=na+nc+1:na+nc+ng;
vd=na+nc+ng+1:na+nc+ng+nd;
Bwo(:,1:2)=Bo(1)*ones(na+nc+ng+nd,2);
Bwo(:,3:4)=Bo(2)*ones(na+nc+ng+nd,2);


Mp=repmat([Mp.*ones(na,1);Mc.*ones(nc,1);Mg.*ones(ng,1)],1,2);
SCw=[Sw(:,1);Cw(:,1);Gw(:,1)];
V0=sum(dVCG.*(1-SCw));
P0=Pi;
Pt0=[P0;Pw(:,1)]';
dtt1=dtt;
while t_flag==1
%for t=1:Ta-1
t=t+1;

if t==1 && dtt==0
    dt=1;
    dtt=1;
else
    dtt=dtt1;
end;
dt1(t+1)=dt;
ft=floor(st);
% if ~isempty(find(CpW(:,t)~=0)) fp=1; else fp=0; end;
fp=1;
    Qf=Qz(:,ft+1);

%     [C,G,A2C,A2G,WonC(:,2),WonG(:,2)]=temp_CG(fll,C2,G2,A2C2,A2G2,WonC1,WonG1,t);
%     [TC,TG,TA2C,TA2G,WonC(:,2),WonG(:,2)]=temp_CG(fll,TC1,TG1,TA2C1,TA2G1,WonC1,WonG1,t);
    
%     [C,TC]=GeoMex(C2,TC1,Pi(na+1:na+nc,1),1);
%     [G,TG]=GeoMex(G2,TG1,Pi(na+nc+1:na+nc+ng,1),1);

     kfw(1:na,1)=Sat_cal(Sw,1,1,as,aw); %water
     kfo(1:na,1)=Sat_cal(Sw,2,1,as,aw); %oil
     
     kfw(na+1:na+nc+ng,1)=Sat_cal([Cw(:,t);Gw(:,t)],1,1,ts,tw); %water
     kfo(na+1:na+nc+ng,1)=Sat_cal([Cw(:,t);Gw(:,t)],2,1,ts,tw); %oil

    [b1gm,b1gc,b1gg,b2gm]=GY_bild(GYData,Pi(1:na,1),Sw(:,1),BZ,na,nc,ng,T_GY,as,aw,mu);
    
    kj=0;
    flag_gim=1;

    while flag_gim==1 && kj<10
        
        kj=kj+1;

        [Clp,Cwp,Cws,A,Bwo,Mp]=SGim(dVCG./Mp(:,1),SCw,Mp,zc,Bwo,Pi,1,P0,va,vc,vg,vd,dt);
        [TL,TW,TP]=Potok_MKT(TM,Pi(1:na,1),kfw(1:na),kfo(1:na),MCp(:,1),as,aw,mu,RC.Arc,mup,fp,kms(1),L,Ke,Ro,A(va));
        [CL,~,~]=Potok_Tube(TC,Pi(vc,1),Cw(:,t),CCp(:,t),PR,mup,fp,kms(2),DATA.Lc,RC.Cr2,RC.Cc2,nc,A(vc));
        [GL,~,~]=Potok_Tube(TG,Pi(vg,1),Gw(:,t),GCp(:,t),PR,mup,fp,kms(3),DATA.Lg,RC.Gr2,RC.Gc2,ng,A(vg));
        
        [A2CL,~,~]=Obmen_T2M(A2C,Pi(1:na,1),Pi(vc,1),Sw(:,1),Cw(:,t),K(:,1),PR,MCp(:,1),CCp(:,t));
        [A2GL,~,~]=Obmen_T2M(A2G,Pi(1:na,1),Pi(vg,1),Sw(:,1),Gw(:,t),K(:,1),PR,MCp(:,1),GCp(:,t));
        
        [W1,W6,W7]=Well_MKT(Wf,Won,Uf(:,ft+1),Sw(:,1),MCp(:,1),aw,as,mu,mup,CpW(:,ft+1),A(va));
        [W1C,W6C,W7C]=Well_MKT(WonC(:,2),WonC(:,1),Uf(WonC(:,3),ft+1),Cw(:,t),CCp(:,t),tw,ts,mu,mup,CpW(WonC(:,3),ft+1),A(vc));
        [W1G,W6G,W7G]=Well_MKT(WonG(:,2),WonG(:,1),Uf(WNG,ft+1),Gw(:,t),GCp(:,t),tw,ts,mu,mup,CpW(WNG,ft+1),A(vg));
        
                
        A1=TL-sparse(Won,Won,W1,na,na)-sparse(1:na,1:na,sum(A2CL,2)+sum(A2GL,2)+Clp(va)+b1gm',na,na);
        C1=CL-sparse(1:nc,1:nc,sum(A2CL,1)+Clp(vc)',nc,nc)-sparse(WonC(:,1),WonC(:,1),W1C,nc,nc);
        G1=GL-sparse(1:ng,1:ng,sum(A2GL,1)+Clp(vg)',ng,ng)-sparse(WonG(:,1),WonG(:,1),W1G,ng,ng);
        
        PwNl=repmat(Pw(:,ft+1),Nl,1);
        b1wm=sparse(Won,ones(1,size(Won,1)),-W1.*PwNl,na,1);
        b1wc=sparse(WonC(:,1),ones(1,size(WonC,1)),-W1C.*Pw(WonC(:,3),ft+1),nc,1);
        b1wg=sparse(WonG(:,1),ones(1,size(WonG,1)),-W1G.*Pw(WNG,ft+1),ng,1);
        
        W2M=sparse(WonM,Won,W1,nw,na);
        W2C=sparse(WonC(:,3),WonC(:,1),W1C,nw,nc);
        W2G=sparse(WNG,WonG(:,1),W1G,nw,ng);
        
        
        b1wm=b1wm.*(sum(W2M(Qf==0,:),1)~=0)';
        b1wc=b1wc.*(sum(W2C(Qf==0,:),1)~=0)';
        b1wg=b1wg.*(sum(W2G(Qf==0,:),1)~=0)';
        %  b1wm
        
        WM1=[W2M,W2C,W2G];
        WM2=WM1';
        W3vec=sparse(WonM,1,W1,nw,1)+sparse(WonC(:,3),1,W1C,nw,1)+sparse(WNG,1,W1G,nw,1);
        WM3=-sparse(1:nw,1:nw,W3vec,nw,nw);
        
        AM=[A1,A2CL,A2GL;A2CL',C1,C2GL;A2GL',C2GL',G1];
        BM=[b1wm;b1wc;b1wg]'+[-b1gm.*GYData.GY_Pz',b1gc,b1gg]-(Clp.*Pi)';
        
        
        WM1=WM1(Qf~=0,:);
        WM2=WM2(:,Qf~=0);
        WM3=WM3(Qf~=0,Qf~=0);
        
        Pt=[BM,Qz(Qf~=0,ft+1)']/[AM,WM2;WM1,WM3];

        flag_gim=sum(abs(Pt(1:na+nc+ng+nd)-Pt0(1:na+nc+ng+nd))./Pt(1:na+nc+ng+nd)>=1e-6)~=0;
        flag_gim=flag_gim*(sum(zc==0)==0);
        Pt0=Pt;
    end
 
    Pi0=Pi;
    Pi(:,1)=Pt(1:na+nc+ng+nd);

    Pw(Qf~=0,ft+1)=Pt(na+nc+ng+1:end);

    %% Водонасыщенность
    SCw=[Sw(:,1);Cw(:,t);Gw(:,t)];
    Sw0=SCw;
    SCp=[MCp(:,1);CCp(:,t);GCp(:,t)];
    
 %[SCw,SCp,NDT(t)]=Sat_fast(SCw,SCp,RC,TC,TG,TA2C,TA2G,TM,Pi(:,1),PR,ndt,Won,Wf,...
 %    Uf(:,t),dt,dVCG,Pw(:,t),WonG,CpW(:,t),WonC,Nl,b2gm,GYData.GY_Pz);
 
 qm(WonM,:)=QBild(W1,W6,W7,Pi(1:na,1),Uf(:,ft+1),Won,dt,Pw(WonM,ft+1));
 qc(WonC(:,3),:)=QBild(W1C,W6C,W7C,Pi(vc,1),Uf(WonC(:,3),ft+1),WonC(:,1),dt,Pw(WonC(:,3),ft+1));
 qg(WonG(:,3),:)=QBild(W1G,W6G,W7G,Pi(vg,1),Uf(WNG,ft+1),WonG(:,1),dt,Pw(WNG,ft+1));
 q=qm+qc+qg;
 Qz1=q(:,1)+q(:,2);

 Qf=Qz1;
 Sw2=Sw;

 [SCw,SCp,NDT(t),Q1,Q2,Qm1,dSS(t)]=Sat_fast_2(SCw,SCp,RC,TC,TG,TA2C,TA2G,Pi(:,1),PR,ndt,Won,Wf,...
     Uf(:,ft+1),dt,dVCG,Pw(:,ft+1),WonG,CpW(:,ft+1),WonC,Nl,CR_rc,CR,Qz1,Qf,Pi0,TW,W6,TP,W7,L,DATA.Lc,DATA.Lg,Ke,Cws,Cwp,TM,A2CL,A2GL,qm);
    
    Sw(:,1)=SCw(1:na);
    Cw(:,t+1)=SCw(na+1:na+nc);
    Gw(:,t+1)=SCw(na+nc+1:end);
    
    MCp(:,1)=SCp(1:na);
    CCp(:,t+1)=SCp(na+1:na+nc);
    GCp(:,t+1)=SCp(na+nc+1:end);

    %% Дебиты
    
    Qm(:,:,t+1)=QBild(W1,W6.*A(Won),W7,Pi(1:na,1),Uf(:,ft+1),Won,dt,Pw(WonM,ft+1));
    Qm(CR_rc.wn,:,t+1)=Qm1;
    Qc(:,:,t+1)=Q1;%QBild(W1C,W6C,W7C,Pi(na+1:na+nc,1),Uf(WonC(:,3),t+1),WonC(:,1),dt,Pw(WonC(:,3),t+1));
    Qg(:,:,t+1)=Q2;%QBild(W1G,W6G,W7G,Pi(na+nc+1:end,1),Uf(WNG,t+1),WonG(:,1),dt,Pw(WNG,t+1));
    PpW(:,t+1)=Pi(Won);
        %% пїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ
%    Ti(:,t+1)=Termal(Ti(:,t),Pi(:,t+1),SCw,Mp,Mc,Mg,NTG,LM,LC,LG,LA2C,LA2G,PR,RC,dVCG,dt,Won,WonM,WonV,WonG,...
 %       TeW(:,t+1),Qm(:,:,t+1),Qc(:,:,t+1),Qg(:,:,t+1),C2GL,TL,CL,GL,TW,CW,GW,A2CL,A2CW,A2GL,A2GW,Wf,NDT(t),t,T,S,BZ);
    
        Ti(:,1)=Ti(:,1);
    if mod(t,1)==0
        j=j+1;
        Pj(:,j)=Pi;
        Swj(:,j)=Sw;
        MCpj(:,j)=MCp;
        Tj(:,j)=Ti;
    end;
    sQo=sum(Qm(:,2,t+1)-Qm(:,3,t+1)+Qm(:,1,t+1))+sum(Qc(:,2,t+1)-Qc(:,3,t+1)+Qc(:,1,t+1))+sum(Qg(:,2,t+1)-Qg(:,3,t+1)+Qg(:,1,t+1));
    %sum(sQo(:))
    dQ(t)=sum(Sw0.*[dV;dVC;dVG])-sum([Sw;Cw(:,t+1);Gw(:,t+1)].*[dV;dVC;dVG])-sum(sQo(:));
    
    dt=vibor_t2(dtt,Pi(1:na),RC,dV,TL,W1,Won,Pw(:,ft+1),na,PR,st,Ta,Sw,Sw2,dt);
    st=st+dt;
    t_flag=st~=Ta;

  
end;

Q=Q2Sut(Qm,Qc,Qg,dt1,Ta);

toc