function [SCw,SCp,ndt,Q1,Q2,Qm,tmp]=Sat_fast_2(SCw,SCp,RC,TC,TG,A2C,A2G,Pi,PR,...
    ndt0,Won,Wf,Uf,dt,dV,Pw,WonG,CpW,WonC,Nl,CR_cr,Qz,Qf,Pi0,TW,W6,TP,W7,L,Lc,Lg,Ke,Cws,Cwp)

na=RC.na;
nc=RC.nc;
ng=RC.ng;

va=1:na;
vc=na+1:na+nc;
vg=na+nc+1:na+nc+ng;

PwNl=repmat(Pw,Nl,1);
 
% aw1=sum(SCw(vc).*dV(vc));

 if isempty(RC.Cr)==0 || isempty(RC.Gr)==0
    [Bc,Bg,SCw(vc),SCw(vg),ndt,Q1,Q2,Qm,dSS]=fun1(RC,Pi,SCw,SCp,PR,TC,TG,A2C,...
        A2G,WonC,WonG,Uf,CpW,Pw,dt,dV,CR_cr,Qz,Qf,ndt0,Pi0,L,Lc,Lg,Ke);
  %  [Bc,Bg,SCw(vc),SCw(vg),ndt,Q1,Q2,Qm,dSS]=fun2(RC,Pi,SCw,SCp,PR,TC,TG,A2C,...
  %      A2G,WonC,WonG,Uf,CpW,Pw,dt,dV,CR_cr,Qz,Qf,ndt0,Pi0,L,Lc,Lg,Ke);
    
    %временно
    Bcp=zeros(size(Bc));
    Bgp=zeros(size(Bc));
    
 else
     Bc=zeros(nc,1);
     Bg=zeros(ng,1);
     
    Bcp=zeros(size(Bc));
    Bgp=zeros(size(Bc));
    
     ndt=1;
     Q1=zeros(size(WonC(:,3),1),5);
     Q2=zeros(size(WonG(:,3),1),5);
     Qm=zeros(0,5);
 end;

%  aw2=sum(SCw(vc).*dV(vc));
%  tmp=-Q1(1)-sum(Bc);
%  AM=TW-sparse(1:na,1:na,sum(TW,2),na,na);%-sparse(Won,Won,W6,na,na);

     b=sparse(Won,ones(1,size(Won,1)),-W6.*(Pi(Won)-PwNl),na,1);
     bp=sparse(Won,ones(1,size(Won,1)),-W7.*(Pi(Won)-PwNl),na,1);

     b(Won(CR_cr.wn))=-(Qm(:,1)+Qm(:,2)-Qm(:,3))/dt;

% sum(Bc)
% sum(Q1)+(Qm(:,1)+Qm(:,2)-Qm(:,3))

 v1=zeros(na,1);
 v1([RC.ACr;RC.AGr])=1;
 r=find(v1==1);
%sum(Bc)

 B=b+sparse(r,ones(sum(v1),1),Bc,na,1)/dt+sparse(r,ones(sum(v1),1),Bg,na,1)/dt-Cwp(va).*(Pi(va)-Pi0(va));
 Bp=bp+sparse(r,ones(sum(v1),1),Bcp,na,1)/dt+sparse(r,ones(sum(v1),1),Bgp,na,1)/dt-SCp(va).*Cwp(va).*(Pi(va)-Pi0(va));

 tmp=sum(B);
% jkghkjh

     SCw_old=SCw(va);
     SCw(va)=SCw(va)+dt*(TW*Pi(va)+B)./Cws(va);

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
