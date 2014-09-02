function [SCw,SCp,ndt]=Sat_fast_CPU(SCw,SCp,RC,TC,TG,A2C,A2G,T,Pi,PR,...
    ndt1,Won,Wf,Uf,dt,dV,Pw,WONG,WNG,CpW)

as=PR.as;
aw=PR.aw;
ts=PR.ts;
tw=PR.tw;
mu=PR.mu;

nc=RC.nc;
nt=RC.nt;
ng=RC.ng;

gPi=gpuArray(Pi);
gSCw=gpuArray(SCw);
gSCp=gpuArray(SCp);
gTC=gpuArray(TC);
gTG=gpuArray(TG);
%A2C=gpuArray(A2C);
%A2G=gpuArray(A2G);
gT=gpuArray(T);

dP=(Pi(RC.rc(:,1))-Pi(RC.rc(:,2)));

i=0; 
j=0;
ndt=1;
b2c=zeros(nt,1);
 
while i==0
  j=j+1;
     [TW,TP]=Potok_MKT_1(gT,gPi(1:nc),gSCw(1:nc),gSCp(1:nc),as,aw,mu,RC.Arc);
     [CW,CP]=Potok_Tube_1(gTC,gPi(nc+1:nc+nt),gSCw(nc+1:nc+nt),gSCp(nc+1:nc+nt),PR,RC.Cr,RC.Cc);
     [GW,GP]=Potok_Tube_1(gTG,gPi(nc+nt+1:end),gSCw(nc+nt+1:end),gSCp(nc+nt+1:end),PR,RC.Gr,RC.Gc);
     
     TW=gather(TW);
     TP=gather(TP);
     CW=gather(CW);
     CP=gather(CP);
     GW=gather(GW);
     GP=gather(GP);

     [A2CW,A2CP]=Obmen_T2M_1(A2C,Pi(1:nc+nt),SCw(1:nc),SCw(nc+1:nc+nt),ts,tw,as,aw,mu,SCp(1:nc),SCp(nc+1:nc+nt),RC.ACr,RC.ACc);
     [A2GW,A2GP]=Obmen_T2M_1(A2G,Pi([1:nc,nc+nt+1:nc+nt+ng]),SCw(1:nc),SCw(nc+nt+1:end),ts,tw,as,aw,mu,SCp(1:nc),SCp(nc+nt+1:end),RC.AGr,RC.AGc);
     
     [W1,W6,W7]=Well_MKT(Wf,Won,Uf,SCw(1:nc),SCp(1:nc),aw,as,mu,CpW);
     [W1G,W6G,W7G]=Well_MKT(WONG(:,1),WONG(:,2),Uf(WNG),SCw(nc+nt+1:end),SCp(nc+nt+1:end),tw,ts,mu,CpW(WNG));
     
     b2=sparse(Won,ones(1,size(Won,1)),W6.*(Pw-Pi(Won)),nc,1);
 %    b2c=zeros(nt,1);
     b2g=sparse(WONG(:,2),ones(1,size(WONG(:,2),1)),W6G.*(Pw(WNG)-Pi(nt+nc+WONG(:,2))),ng,1);

     b2p=sparse(Won,ones(1,size(Won,1)),W7.*(Pw-Pi(Won)),nc,1);
     b2gp=sparse(WONG(:,2),ones(1,size(WONG(:,2),1)),W7G.*(Pw(WNG)-Pi(nt+nc+WONG(:,2))),ng,1);

     V=[TW;CW;A2CW;A2CW;GW;A2GW;A2GW].*dP;
     V1=accumarray(RC.rc(:,2),V);
     
     Vp=[TP;CP;A2CP;A2CP;GP;A2GP;A2GP].*dP;
     V1p=accumarray(RC.rc(:,2),Vp);

    % Vn=sparse(rr,cc,V,nc+nt,nc+nt);
    % V1=sum(Vn,1)'-vb2;
     SCw_old=SCw;
%         V1=gpuArray(V1);
%         V1p=gpuArray(V1p);
        
     if j==1
         SCw1=SCw+dt*(V1+[b2;b2c;b2g])./dV;%SCw+
         SCp1=(SCw_old.*SCp+dt*(V1p+[b2p;b2c;b2gp])./dV);%SCw+
         SCp1(SCw1~=0)=SCp1(SCw1~=0)./SCw1(SCw1~=0);
         
         mr=sum(SCw1>1)+sum(SCw1<0)+sum(SCp1>1)+sum(SCp1<0);
         
         if mr~=0
             ndt=ceil(max(abs(SCw1)))*ndt1;
             SCw=SCw+dt/ndt*(V1+[b2;b2c;b2g])./dV;
             SCp=(SCw_old.*SCp+dt/ndt*(V1p+[b2p;b2c;b2gp])./dV);
             SCp(SCw~=0)=SCp(SCw~=0)./SCw(SCw~=0);
         else
             SCw=SCw1;
             SCp=SCp1;
         end;
     else
         SCw=SCw+dt/ndt*(V1+[b2;b2c;b2g])./dV;
         SCp=(SCw_old.*SCp+dt/ndt*(V1p+[b2p;b2c;b2gp])./dV);
         SCp(SCw~=0)=SCp(SCw~=0)./SCw(SCw~=0);
     end;
    
     SCw=SCw.*(SCw>=0).*(SCw<=1)+(SCw>1);
     SCp=SCp.*(SCp>=0).*(SCp<=1)+(SCp>1);
     
  i=j>=ndt; 

end;

% SCw=gather(SCw);
% SCp=gather(SCp);