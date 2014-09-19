function [XY,K,Z,Pi,Sw,MCp,p,Q]=Sim_MKT(PR,C,A2C,G,A2G,dVC,dVG,DATA,WData)
tic

K=DATA.gK;
Mp=DATA.gMp;
P=DATA.gP;
Sw=DATA.gSw;
MCp=DATA.gCp;
XY=DATA.XY;
H=DATA.gH;
Z=DATA.gZ;
Won=DATA.Won;
WnoG=DATA.WonG;
%WNG=DATA.WNG;
WNG=DATA.WonG(:,3);
WonV=DATA.WonV;

Pw=WData.Pw;
Uf=WData.Uf;
CpW=WData.CpW;

Nl=PR.Nl;
as=PR.as;
aw=PR.aw;
ts=PR.ts;
tw=PR.tw;
mu=PR.mu;
Ta=PR.Ta;
dt=PR.dt;
ndt=PR.ndt;

%Q=zeros(size(Uf,1),size(Uf,2),4);
Qc=zeros(size(WonV,1),size(Uf,2),5);
Qg=zeros(size(WNG,1),size(Uf,2),5);
K(:)=mean(K(:));

[A]=MR_Prop(XY,Nl);

[L,B,S,H1]=Geome3(A,XY,Z,H);
Wf=KWell(K,H,S,Won);

p=symrcm(A);
A=A(p,p);  L=L(p,p);
S=S(p,p);  B=B(p,p); 
H1=H1(p,p);
[r,c]=find(L);

XY=repmat(XY,Nl,1);
K=K(:); Mp=Mp(:); Sw=Sw(:); H=H(:); Z=Z(:); P=P(:); MCp=MCp(:);
K=K(p); Mp=Mp(p); Sw=Sw(p); XY=XY(p,:); H=H(p); Z=Z(p); P=P(p); MCp=MCp(p);

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

nc=size(XY,1);  nt=size(C,1);  ng=size(G,1);
C2G=sparse(nt,ng);
C2GL=C2G;

figure(98),subplot(2,4,7),spy([A,A2C,A2G;A2C',C,C2G;A2G',C2G',G]);

Pi(1:nc,1)=P;
Pi(nc+1:nc+nt,1)=P(1:nt);
Pi(nc+nt+1:nc+nt+ng,1)=P(1:ng);

[Ke,dV]=KH2Mat(K,H,Mp,S,r,c); 
b1_2_C=zeros(nt,1);

Cw(:,1)=Sw(1:nt,1);
Gw(:,1)=Sw(1:ng,1);

CCp(:,1)=zeros(nt,1);
GCp(:,1)=zeros(ng,1);

dVCG=[dV;dVC;dVG];

[T,TC,TG,TA2C,TA2G,RC]=Pre_fast(A,C,A2C,A2G,G,Ke,L,B,H1,K);

for t=1:Ta
    [TL,TW]=Potok_MKT(T,Pi(1:nc,t),Sw(:,t),MCp(:,t),as,aw,mu,RC.Arc);
    [CL,CW]=Potok_Tube(C,Pi(nc+1:nc+nt,t),Cw(:,t),CCp(:,t),PR);
    [GL,GW]=Potok_Tube(G,Pi(nc+nt+1:end,t),Gw(:,t),GCp(:,t),PR);

    [A2CL,A2CW]=Obmen_T2M(A2C,Pi(:,t),Sw(:,t),Cw(:,t),K,PR,MCp(:,t),CCp(:,t));
    [A2GL,A2GW]=Obmen_T2M(A2G,Pi(:,t),Sw(:,t),Gw(:,t),K,PR,MCp(:,t),GCp(:,t));

    %[gTL,gTW]=Potok_MKT_GPU(L,B,Ke,Pi(:,t),Sw(:,t),H1,r,c,as,aw);
    
    [W1,W6,W7]=Well_MKT(Wf,Won,Uf(:,t),Sw(:,t),MCp(:,t),aw,as,mu,CpW(:,t));
    [W1C,W6C,W7C]=Well_MKT(WonV(:,2),WonV(:,1),Uf(WonV(:,3),t),Cw(:,t),CCp(:,t),tw,ts,mu,CpW(WonV(:,3),t));
    [W1G,W6G,W7G]=Well_MKT(WnoG(:,1),WnoG(:,2),Uf(WNG,t),Gw(:,t),GCp(:,t),tw,ts,mu,CpW(WNG,t));
    
    A1=TL-sparse(Won,Won,W1,nc,nc)-sparse(1:nc,1:nc,sum(A2CL,2),nc,nc)-sparse(1:nc,1:nc,sum(A2GL,2),nc,nc);
    C1=CL-sparse(1:nt,1:nt,sum(A2CL,1),nt,nt)-sparse(WonV(:,1),WonV(:,1),W1C,nt,nt);
    G1=GL-sparse(1:ng,1:ng,sum(A2GL,1),ng,ng)-sparse(WnoG(:,2),WnoG(:,2),W1G,ng,ng);

    b1=sparse(Won,ones(1,size(Won,1)),-W1.*Pw(:,t),nc,1);

%     size(-W1G)
%     size(Pw(WonV(:,3),t))
% nt

    b1_2_C=sparse(WonV(:,1),ones(1,size(WonV,1)),-W1C.*Pw(WonV(:,3),t),nt,1);
    b1_2_G=sparse(WnoG(:,2),ones(1,size(WnoG,1)),-W1G.*Pw(WNG,t),ng,1);

    Pi(:,t+1)=[b1;b1_2_C;b1_2_G]'/[A1,A2CL,A2GL;A2CL',C1,C2GL;A2GL',C2GL',G1];
    
    
    b2=sparse(Won,ones(1,size(Won,1)),W6.*Pw(:,t),nc,1);
    b2_2_C=sparse(WonV(:,1),ones(1,size(WonV(:,1),1)),W6C.*Pw(WonV(:,3),t),nt,1);
    b2_2_G=sparse(WnoG(:,2),ones(1,size(WnoG(:,2),1)),W6G.*Pw(WNG,t),ng,1);
    
    b2p=sparse(Won,ones(1,size(Won,1)),W7.*Pw(:,t),nc,1);
    b2p_2_C=sparse(WonV(:,1),ones(1,size(WonV(:,1),1)),W7C.*Pw(WonV(:,3),t),nt,1);
    b2p_2_G=sparse(WnoG(:,2),ones(1,size(WnoG(:,2),1)),W7G.*Pw(WNG,t),ng,1);
    
    SCw=[Sw(:,t);Cw(:,t);Gw(:,t)];
    SCp=[MCp(:,t);CCp(:,t);GCp(:,t)];

 [SCw,SCp,NDT(t)]=Sat_fast(SCw,SCp,RC,TC,TG,TA2C,TA2G,T,Pi(:,t+1),PR,ndt,Won,Wf,...
     Uf(:,t),dt,dVCG,Pw(:,t),WnoG,CpW(:,t),WonV);
    
%       [SCw,SCp,NDT(t)]=Sat_fast_CPU(SCw,SCp,RC,TC,TG,TA2C,TA2G,T,Pi(:,t),PR,ndt,Won,Wf,...
%           Uf(:,t),dt,dVCG,Pw(:,t),WONG,WNG,CpW(:,t));
    
    Sw(:,t+1)=SCw(1:nc);
    Cw(:,t+1)=SCw(nc+1:nc+nt);
    Gw(:,t+1)=SCw(nc+nt+1:end);
    
    MCp(:,t+1)=SCp(1:nc);
    CCp(:,t+1)=SCp(nc+1:nc+nt);
    GCp(:,t+1)=SCp(nc+nt+1:end);
   
    Sw(:,t+1)=Sw(:,t+1).*(Sw(:,t+1)>=0).*(Sw(:,t+1)<=1)+(Sw(:,t+1)>1);
    Cw(:,t+1)=Cw(:,t+1).*(Cw(:,t+1)>=0).*(Cw(:,t+1)<=1)+(Cw(:,t+1)>1);
    Gw(:,t+1)=Gw(:,t+1).*(Gw(:,t+1)>=0).*(Gw(:,t+1)<=1)+(Gw(:,t+1)>1);
    
    Qm(:,t,:)=QBild(W1,W6,W7,b1,b2,b2p,Pi(1:nc,t+1),Uf(:,t),Won,nc,dt);
    Qc(:,t,:)=QBild(W1C,W6C,W7C,b1_2_C,b2_2_C,b2p_2_C,Pi(nc+1:nc+nt,t+1),Uf(WonV(:,3),t),WonV(:,1),nt,dt);
    Qg(:,t,:)=QBild(W1G,W6G,W7G,b1_2_G,b2_2_G,b2p_2_G,Pi(nc+nt+1:end,t+1),Uf(WNG,t),WnoG(:,2),ng,dt);

end;
% MCp(MCp(:,end)~=0,:)
 %Sw(Sw(:,end)~=1,:)
% [r,c]=find(MCp(:,end)~=0)
% CCp

Q(:,:,1)=[Qm(:,:,1);Qc(:,:,1);Qg(:,:,1)];
Q(:,:,2)=[Qm(:,:,2);Qc(:,:,2);Qg(:,:,2)];
Q(:,:,3)=[Qm(:,:,3);Qc(:,:,3);Qg(:,:,3)];
Q(:,:,4)=[Qm(:,:,4);Qc(:,:,4);Qg(:,:,4)];

% max(SDF,[],1)
% max(SDF1,[],1)
%sum(SDF(:)>1)
%sum(SDF1(:)>1)

sum(NDT)-Ta
toc