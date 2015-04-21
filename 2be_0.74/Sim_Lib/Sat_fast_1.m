function [Sw,Cp,ndt,Q1,Q2,Qm,tmp]=Sat_fast_1(Sw,Cp,RC,TC,TG,A2C,A2G,Pi,PR,...
    ndt0,Won,Uf,dt,dV,Pw,WonG,CpW,WonC,Nl,CR_ind,Qz,Qf,Pi0,TL,W1,TW,W6,TP,...
    W7,L,Lc,Lg,Ke,Cws,Cwp,BLGY_GIM,Qzm1,WonM,nw,b1gm,Clp)

na=RC.na;
nc=RC.nc;
ng=RC.ng;

va=1:na;
vc=na+1:na+nc;
vg=na+nc+1:na+nc+ng;

NwC=size(WonC,1);

PwNl=repmat(Pw,Nl,1);
     v1=zeros(na,1);
     v1([RC.ACr;RC.AGr])=1;
     r=find(v1==1);
 
% aw1=sum(SCw(vc).*dV(vc));

 if isempty(RC.Cr)==0 || isempty(RC.Gr)==0
  bwc=zeros(na,size(CR_ind,1));
  bwg=zeros(na,size(CR_ind,1));
  blc=zeros(na,size(CR_ind,1));
  blg=zeros(na,size(CR_ind,1));
  Qm1=zeros(size(Qz,1));
  
     for i=1:size(CR_ind,1)
      CR_cr=CR_ind{i,1};
      TC=CR_ind{i,2};
      TG=CR_ind{i,3};
      A2C=CR_ind{i,4};
      A2G=CR_ind{i,5};
      RC_p=CR_ind{i,6};
      vc_p=CR_ind{i,7};
      rv=find(CR_cr.v2==1);
      
      RC.Cr2=RC_p.Cr2;
      RC.Cc2=RC_p.Cc2;
      RC.Gr2=RC_p.Gr2;
      RC.Gc2=RC_p.Gc2;
      RC.ACc=RC_p.ACc;
      RC.AGc=RC_p.AGc;
      RC.ACr=RC_p.ACr;
      RC.AGr=RC_p.AGr;
      RC.nc=size(vc_p,2);
      wonC=CR_ind{i,8};
      WonG=CR_ind{i,9};            


     [bwc_1,bwg_1,blc_1,blg_1,Sw(vc_p),Sw(vg),ndt(i),Q1_i,Q2,Qm_i,dSS]=fun1(RC,Pi([va,vc_p]),Sw([va,vc_p]),Cp([va,vc_p]),PR,TC,TG,A2C,...
         A2G,wonC,WonG,Uf,CpW,Pw,dt,dV,CR_cr,Qz,Qf,ndt0,Pi0([va,vc_p]),L,Lc,Lg,Ke);

      bwc(:,i)=sparse(rv,1,bwc_1,na,1);
    %  bwg(:,i)=sparse(r,1,bwg_1,na,1);
      blc(:,i)=sparse(rv,1,blc_1,na,1);
     % blg(:,i)=sparse(r,1,blg_1,na,1);
%     Qm_i
      Qm1(CR_cr.wn,:)=Qm_i;
%      i
      QW1(wonC(:,3),:)=Q1_i;
     end;

     for i=1:NwC 
        Q1(i,:)=QW1(WonC(i,3),:);
     end
     Blc=sum(blc,2);
     Blg=sum(blg,2);
     Qm=Qm1;
     
     
     WM1=sparse(WonM,Won,W1,nw,na);
     WM2=WM1';
     W3vec=sparse(WonM,1,W1,nw,1);
     WM3=-sparse(1:nw,1:nw,W3vec,nw,nw);
  
     Qf=Qzm1<0;
     
     b1wm=sparse(Won,ones(1,size(Won,1)),-W1.*PwNl,na,1);
     b1wm=b1wm.*(sum(WM1(Qf==0,:),1)~=0)';
     bl=b1wm'+BLGY_GIM;
     
     Bl=bl'-Blc/dt;%-sparse(r,ones(sum(v1),1),Blg,na,1)/dt;
     
     A1=TL-sparse(Won,Won,W1,na,na)-sparse(1:na,1:na,Clp(va)+b1gm',na,na);
     
     WM1=WM1(Qf~=0,:);
     WM2=WM2(:,Qf~=0);
     WM3=WM3(Qf~=0,Qf~=0);
     Qzm1(CR_cr.wn)=(Qm_i(:,1)+Qm_i(:,2))/dt;
     Pt=[Bl',Qzm1(Qf~=0)']/[A1,WM2;WM1,WM3];
     
     Pi(va)=Pt(va)';
     PwNl(Qf~=0)=Pt(na+1:end);
     
     %временно
     Bpc=zeros(size(Blc));
     Bpg=zeros(size(Blc));
    
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

Bwc=sum(bwc,2);
Bwg=sum(bwg,2);


     bw=sparse(Won,ones(1,size(Won,1)),-W6.*(Pi(Won)-PwNl),na,1);
     bp=sparse(Won,ones(1,size(Won,1)),-W7.*(Pi(Won)-PwNl),na,1);
     
     Bw=bw+Bwc/dt-Cwp(va).*(Pi(va)-Pi0(va));%+sparse(r,ones(sum(v1),1),Bwg,na,1)/dt
     Bp=bp+Bpc/dt-Cp(va).*Cwp(va).*(Pi(va)-Pi0(va)); %+sparse(r,ones(sum(v1),1),Bpg,na,1)/dt
     tmp=sum(Bw);    

     Sw_old=Sw(va);
     Sw(va)=Sw(va)+dt*(TW*Pi(va)+Bw)./Cws(va);

     Cp(va)=Sw_old.*Cp(va)+dt*(TP*Pi(va)+Bp)./Cws(va);
     v0=Sw(va)==0;
     Cp(va(v0==0))=Cp(va(v0==0))./Sw(va(v0==0));

     Sw=Sw.*(Sw>=0).*(Sw<=1)+(Sw>1);
     Cp=Cp.*(Cp>=0).*(Cp<=1)+(Cp>1);