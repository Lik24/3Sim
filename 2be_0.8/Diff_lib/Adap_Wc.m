function [dJ1,Dd]=Adap_Wc(PR,C,A2C,G,A2G,dVC,dVG,DATA,WData,GYData,fll,Q,Qf,Pp,Swp,a,dJp,Dp,it)

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

Nl=PR.Nl;  as=PR.as; aw=PR.aw; ts=PR.ts; tw=PR.tw; mu=PR.mu;
Ta=PR.Ta;  dt=PR.dt; ndt=PR.ndt; zc=PR.zc; Bo=PR.Bo; 


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

%Wf=KWell_2(KX,H,S,L,B,Won,r,c,WData.Doly,WData.r0);
Wf=KWell_4(KX,H,S,L,B,Won,r,c,WData.Doly,WData.r0,XY);

% Wf(1)=Wf(9);
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

na=size(XY,1);  nc=size(C,1);  ng=size(G,1);    nw=size(Qz,1);
C2G=sparse(nc,ng);
C2GL=C2G;
WonM=repmat(1:nw,1,Nl);

[Ke,Ke_gy,dV]=KH2Mat(K,H,Mp,S,r,c,rz,cz,GYData.GY_Kz,BZ); 

Mc=ones(nc,1);
Mg=ones(ng,1);

CCp(:,1)=zeros(nc,1);
GCp(:,1)=zeros(ng,1);

dVCG=[dV;dVC;dVG];

[TM,TC,TG,TA2C,TA2G,RC,T_GY]=Pre_fast(A,C,A2C,A2G,G,Ke,L,B,H1,K(:,1),Ke_gy,BZ,rz,cz);
Cw(:,1)=Sw(RC.ACr,1);
Gw(:,1)=Sw(RC.AGr,1);
Pi(1:na,1)=P;
Pi(na+1:na+nc,1)=P(RC.ACr);
Pi(na+nc+1:na+nc+ng,1)=P(RC.AGr);

Ti(1:na,1)=T;
Ti(na+1:na+nc,1)=T(RC.ACr);
Ti(na+nc+1:na+nc+ng,1)=T(RC.AGr);
Ksi1=zeros(na,Ta+1); Ksi2=zeros(na,Ta+1);
%[LM,LC,LG,LA2C,LA2G]=Pre_fast(A,C,A2C,A2G,G,1,L,B,H1,ones(size(K(:,1))));
j=0;
PpW=zeros(size(Uf));

[Img2,CR_rc]=Pre_Crack(RC,na,TM,A2C,A2G,Wf,Won,WonM,nw);


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
D=dV/dt;
D1=diag(D);
% t=5;
% [dW1,dW6,W7]=Well_MKT(Wf,Won,Uf(:,t),Swp(1:na,t-1),MCp(:,1),aw,as,mu,CpW(:,t),2);  %производные
% [dTLp,dTWp] = dT_dx2(TM,Pp(1:na,t-1:t),Swp(1:na,t-1),MCp(:,1),as,aw,mu,RC.Arc); %производные
% Kfw=Sat_cal(Swp(1:na,t-1),1,2,as,aw); %water
% Kfo=Sat_cal(Swp(1:na,t-1),2,2,as,aw); %oil
for t=Ta:-1:2
  
[TL,TW]=Potok_MKT(TM,Pp(1:na,t-1),Swp(1:na,t-1),MCp(:,1),as,aw,mu,RC.Arc,1); 
[W1,W6,W7]=Well_MKT(Wf,Won,Uf(:,t),Swp(1:na,t-1),MCp(:,1),aw,as,mu,CpW(:,t),1);

[dTLp,dTWp] = dT_dx2(TM,Pp(1:na,t-1:t),Swp(1:na,t-1),MCp(:,1),as,aw,mu,RC.Arc); %производные
[dW1,dW6,W7]=Well_MKT(Wf,Won,Uf(:,t),Swp(1:na,t-1),MCp(:,1),aw,as,mu,CpW(:,t),2);  %производные

% dF1=dTLp+sparse(Won,Won,dW1.*(Pw(:,t)-Pp(Won,t)),na,na);
% dF2=D1+dTWp+sparse(Won,Won,dW6.*(Pw(:,t)-Pp(Won,t)),na,na);
[dF1 ] = dF_dx2(dTLp,dW1,Pw(:,t),Pp(Won,t),Won,na,0,1 );
[dF2 ] = dF_dx2(dTWp,dW6,Pw(:,t),Pp(Won,t),Won,na,D1,2 );

f=Q(:,t)-Qf(:,t);

A1=TL-sparse(Won,Won,W1,na,na);
A2=TW-sparse(Won,Won,W6,na,na);
%Решение сопряженной задачи
Ksi2(:,t)=((dF1)'*Ksi1(:,t+1)+(dF2)'*Ksi2(:,t+1))./D;
Ksi1(:,t)=A1'\(-A2'*Ksi2(:,t)-(2)*sparse(Won,WonM,W1,na,nw)*f);
%Нахождение градиента целевого функционала
%  [dF1u] = -dF_du(W1,a,Pp(:,t),sparse(Won,1,Pw(:,t),na,1),1);
%  [dF2u] = -dF_du(W6,a,Pp(:,t),sparse(Won,1,Pw(:,t),na,1),1);
%  
%  dF1u=dF1u.*sparse(Won,WonM,ones(nw,1),na,nw);
%  dF2u=dF2u.*sparse(Won,WonM,ones(nw,1),na,nw);
%  [dfu] = dF_du(W1,a,Pp(:,t),sparse(Won,1,Pw(:,t),na,1),2);
[dF1u] = -dF_du(W1,a,Pp(Won,t),Pw(:,t),1,Won,WonM,na,nw);
[dF2u] = -dF_du(W6,a,Pp(Won,t),Pw(:,t),1,Won,WonM,na,nw);
[dfu] = dF_du(W1,a,Pp(Won,t),Pw(:,t),2,Won,WonM,na,nw);
 dJ(:,t)=-((dF1u)'*Ksi1(:,t)+(dF2u)'*Ksi2(:,t)+(2)*(dfu)'*f)*dt;
% if t==710
%     sadas
% end
end
dJ1=sum(dJ,2);
if it>1
    y=(sum(dJ1.*dJ1)/sum(dJp.*dJp));
    Dd=dJ1+y*Dp;
else
    Dd=dJ1;
    y=0;
end
%b=mean(abs(dJ1(find(dJ1~=0))));
%toc;
end
