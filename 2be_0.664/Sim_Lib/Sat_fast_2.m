function [SCw,SCp,ndt,Q1,Q2,Qm,tmp]=Sat_fast_2(SCw,SCp,RC,TC,TG,A2C,A2G,Pi,PR,...
    ndt0,Won,Wf,Uf,dt,dV,Pw,WonG,CpW,WonC,Nl,CR_cr,CR,Qz,Qf,Pi0,TL,W1,TW,W6,TP,...
    W7,L,Lc,Lg,Ke,Cws,Cwp,TM,A2CL,A2GL,Qm0,BLGY_GIM,Qzm1,WonM,nw)

na=RC.na;
nc=RC.nc;
ng=RC.ng;

va=1:na;
vc=na+1:na+nc;
vg=na+nc+1:na+nc+ng;

PwNl=repmat(Pw,Nl,1);
 
% aw1=sum(SCw(vc).*dV(vc));

 if isempty(RC.Cr)==0 || isempty(RC.Gr)==0
    [Bwc,Bwg,Blc,Blg,SCw(vc),SCw(vg),ndt,Q1,Q2,Qm,dSS]=fun1(RC,Pi,SCw,SCp,PR,TC,TG,A2C,...
        A2G,WonC,WonG,Uf,CpW,Pw,dt,dV,CR_cr,Qz,Qf,ndt0,Pi0,L,Lc,Lg,Ke);
  %  [Bc,Bg,SCw(vc),SCw(vg),ndt,Q1,Q2,Qm,dSS]=fun2(RC,Pi,SCw,SCp,PR,TC,TG,A2C,...
   %     A2G,WonC,WonG,Uf,CpW,Pw,dt,dV,CR_cr,Qz,Qf,ndt0,Pi0,L,Lc,Lg,Ke,CR,TM);
    
Blc0=(A2CL*Pi(na+1:na+nc)-sum(A2CL,2).*Pi(1:na))*dt;
Blg0=(A2GL*Pi(na+nc+1:na+nc+ng)-sum(A2GL,2).*Pi(1:na))*dt;
    %временно
    Bpc=zeros(size(Bwc));
    Bpg=zeros(size(Bwc));
    
 else
     Bwc=zeros(nc,1);
     Bwg=zeros(ng,1);
     
    Bpc=zeros(size(Bwc));
    Bpg=zeros(size(Bwc));
    
     ndt=1;
     Q1=zeros(size(WonC(:,3),1),5);
     Q2=zeros(size(WonG(:,3),1),5);
     Qm=zeros(0,5);
 end;


%  aw2=sum(SCw(vc).*dV(vc));
%  tmp=-Q1(1)-sum(Bc);
%  AM=TW-sparse(1:na,1:na,sum(TW,2),na,na);%-sparse(Won,Won,W6,na,na);
WM1=sparse(WonM,Won,W1,nw,na);
WM2=WM1';
W3vec=sparse(WonM,1,W1,nw,1);
WM3=-sparse(1:nw,1:nw,W3vec,nw,nw);

Qf=Qzm1<0;

     b1wm=sparse(Won,ones(1,size(Won,1)),-W1.*PwNl,na,1);
     b1wm=b1wm.*(sum(WM1(Qf==0,:),1)~=0)';
     BM=b1wm'+BLGY_GIM;
     bl=BM';


%     vgmnoj=Qm(:,1)/Qm0(CR_cr.wn,1);
%Qm=Qm./vgmnoj;

     bw(Won(CR_cr.wn))=-(Qm(:,1)+Qm(:,2)-Qm(:,3))/dt;

% sum(Bc)
% sum(Q1)+(Qm(:,1)+Qm(:,2)-Qm(:,3))

 v1=zeros(na,1);
 v1([RC.ACr;RC.AGr])=1;
 r=find(v1==1);
%[Blc-Blc0(r)]
%vmoj=Blc./Blc0(r);
%Bc=Bc./vmoj;

Bl=bl-sparse(r,ones(sum(v1),1),Blc,na,1)/dt-sparse(r,ones(sum(v1),1),Blg,na,1)/dt;

A1=TL-sparse(Won,Won,W1,na,na);
         
       
WM1=WM1(Qf~=0,:);
WM2=WM2(:,Qf~=0);
WM3=WM3(Qf~=0,Qf~=0);
Qzm1(CR_cr.wn)=Qm(:,1)+Qm(:,2);
Pt=[Bl',Qzm1(Qf~=0)'/dt]/[A1,WM2;WM1,WM3];

Pi(va)=Pt(va)';
PwNl(Qf~=0)=Pt(na+1:end);

     bw=sparse(Won,ones(1,size(Won,1)),-W6.*(Pi(Won)-PwNl),na,1);
     bp=sparse(Won,ones(1,size(Won,1)),-W7.*(Pi(Won)-PwNl),na,1);
     
     Bw=bw+sparse(r,ones(sum(v1),1),Bwc,na,1)/dt+sparse(r,ones(sum(v1),1),Bwg,na,1)/dt-Cwp(va).*(Pi(va)-Pi0(va));
     Bp=bp+sparse(r,ones(sum(v1),1),Bpc,na,1)/dt+sparse(r,ones(sum(v1),1),Bpg,na,1)/dt-SCp(va).*Cwp(va).*(Pi(va)-Pi0(va));
 tmp=sum(Bw);    
SD1=[TW*Pi(va)+Bw];%,full(sparse(r,ones(sum(v1),1),Bc,na,1)/dt),full(b)]
%SD1(598,:)
%r

%jhjh
     SCw_old=SCw(va);
     SCw(va)=SCw(va)+dt*(TW*Pi(va)+Bw)./Cws(va);

     SCp(va)=SCw_old.*SCp(va)+dt*(TP*Pi(va)+Bp)./Cws(va);
     v0=SCw(va)==0;
     SCp(va(v0==0))=SCp(va(v0==0))./SCw(va(v0==0));

     SCw=SCw.*(SCw>=0).*(SCw<=1)+(SCw>1);
     SCp=SCp.*(SCp>=0).*(SCp<=1)+(SCp>1);
    % ndt=1;
%      SCw'
%      hgj
%Bl=sum(Bcl);
%Bl
