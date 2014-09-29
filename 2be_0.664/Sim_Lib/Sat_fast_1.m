function [SCw,SCp,ndt]=Sat_fast_1(SCw,SCp,RC,TC,TG,A2C,A2G,T,Pi,PR,...
    ndt1,Won,Wf,Uf,dt,dV,Pw,WonG,CpW,WonV,Nl,b2w,P_gy)

as=PR.as;
aw=PR.aw;
ts=PR.ts;
tw=PR.tw;
mu=PR.mu;

na=RC.na;
nc=RC.nc;
ng=RC.ng;

dP=(Pi(RC.rc(:,1))-Pi(RC.rc(:,2)));

i=0; 
j=0;
ndt=1;
%b2c=zeros(nt,1);

[vPa,vPc1,vPc2,vPg1,vPg2]=pre_potok(Pi,RC);
iL=ones(na,1);
sKfw=[];

while i==0
  j=j+1;
  
     [TW1,TP1,mr,sKfw]=Potok_Vib(SCw(1:na),SCp(1:na),RC,Pi,as,aw,iL,T,mu,sKfw);
     
     [TW,TP]=Potok_MKT_1(T,vPa,SCw(1:na),SCp(1:na),as,aw,mu,RC.Arc);
%sum(TW1-TW(mr==1))
    % [TW1,TP1]=Potok_MKT_1(T(AV==1),vPa(AV==1,:),SWM(iL==1),CPM(iL==1),as,aw,mu,[mr,mc]);

%      size(TW1(:,1))/size(TW(:,1))
% 
%      if j==50
%          ghj
%      end;

     [CW,CP]=Potok_Tube_1(TC,Pi(na+1:na+nc),vPc1,vPc2,SCw(na+1:na+nc),SCp(na+1:na+nc),PR,RC.Cr,RC.Cc);
     [GW,GP]=Potok_Tube_1(TG,Pi(na+nc+1:end),vPg1,vPg2,SCw(na+nc+1:end),SCp(na+nc+1:end),PR,RC.Gr,RC.Gc);

     [A2CW,A2CP]=Obmen_T2M_1(A2C,Pi(1:na+nc),SCw(1:na),SCw(na+1:na+nc),ts,tw,as,aw,mu,SCp(1:na),SCp(na+1:na+nc),RC.ACr,RC.ACc);
     [A2GW,A2GP]=Obmen_T2M_1(A2G,Pi([1:na,na+nc+1:na+nc+ng]),SCw(1:na),SCw(na+nc+1:end),ts,tw,as,aw,mu,SCp(1:na),SCp(na+nc+1:end),RC.AGr,RC.AGc);

     % full([sum(A2CW(:)),sum(A2GW(:))])
      
     [W1,W6,W7]=Well_MKT(Wf,Won,Uf,SCw(1:na),SCp(1:na),aw,as,mu,CpW);
     [W1C,W6C,W7C]=Well_MKT(WonV(:,2),WonV(:,1),Uf(WonV(:,3)),SCw(na+1:na+nc),SCp(na+1:na+nc),tw,ts,mu,CpW(WonV(:,3)));
     [W1G,W6G,W7G]=Well_MKT(WonG(:,2),WonG(:,1),Uf(WonG(:,3)),SCw(na+nc+1:end),SCp(na+nc+1:end),tw,ts,mu,CpW(WonG(:,3)));
     
     PwNl=repmat(Pw,Nl,1);
     
     b=sparse(Won,ones(1,size(Won,1)),W6.*(PwNl-Pi(Won)),na,1)+b2w'.*(P_gy-Pi(1:na));
     bc=sparse(WonV(:,1),ones(1,size(WonV,1)),W6C.*(Pw(WonV(:,3))-Pi(na+WonV(:,1))),nc,1);
     bg=sparse(WonG(:,1),ones(1,size(WonG,1)),W6G.*(Pw(WonG(:,3))-Pi(nc+na+WonG(:,1))),ng,1);

     bp=sparse(Won,ones(1,size(Won,1)),W7.*(PwNl-Pi(Won)),na,1);
     bcp=sparse(WonV(:,1),ones(1,size(WonV(:,1),1)),W7C.*(Pw(WonV(:,3))-Pi(na+WonV(:,1))),nc,1);
     bgp=sparse(WonG(:,1),ones(1,size(WonG(:,1),1)),W7G.*(Pw(WonG(:,3))-Pi(nc+na+WonG(:,1))),ng,1);

     V=[TW;CW;A2CW;A2CW;GW;A2GW;A2GW].*dP;
     V1=accumarray(RC.rc(:,2),V,[na+nc+ng,1]);
     
     dP1=dP(1:size(T,1));
     dP1=dP1(mr==1);
     dP2=dP(size(T,1)+1:end);
     dP3=[dP1;dP2];
     
     Rrc=RC.rc(1:size(T,1),2);
     Rrc=Rrc(mr==1);
     Rrc2=RC.rc(size(T,1)+1:end,2);
     Rrc3=[Rrc;Rrc2];
     
     Vv=[TW1;CW;A2CW;A2CW;GW;A2GW;A2GW].*dP3;
%      size(Rrc3)
%      size(Vv)
     V2=accumarray(Rrc3,Vv);
     
   
     Vp=[TP;CP;A2CP;A2CP;GP;A2GP;A2GP].*dP;
     V1p=accumarray(RC.rc(:,2),Vp,[na+nc+ng,1]);


     %V2=accumarray(RC.rc(:z,2),V,[na+nc+ng,1]);
    % Vn=sparse(rr,cc,V,nc+nt,nc+nt);
    % V1=sum(Vn,1)'-vb2;
     SCw_old=SCw;

        
     if j==1
         SCw1=SCw+dt*(V1+[b;bc;bg])./dV;%SCw+
         SCp1=(SCw_old.*SCp+dt*(V1p+[bp;bcp;bgp])./dV);%SCw+
         SCp1(SCw1~=0)=SCp1(SCw1~=0)./SCw1(SCw1~=0);
                  

         mr=sum(SCw1>1)+sum(SCw1<0)+sum(SCp1>1)+sum(SCp1<0);
         
         if mr~=0
             ndt=ceil(max(abs(SCw1)))*ndt1;
             SCw=SCw+dt/ndt*(V1+[b;bc;bg])./dV;
             SCp=(SCw_old.*SCp+dt/ndt*(V1p+[bp;bcp;bgp])./dV);
             SCp(SCw~=0)=SCp(SCw~=0)./SCw(SCw~=0);
         else
             SCw=SCw1;
             SCp=SCp1;
         end;
     else
         SCw=SCw+dt/ndt*(V1+[b;bc;bg])./dV;
         SCp=(SCw_old.*SCp+dt/ndt*(V1p+[bp;bcp;bgp])./dV);
         SCp(SCw~=0)=SCp(SCw~=0)./SCw(SCw~=0);
     end;

        
     SCw=SCw.*(SCw>=0).*(SCw<=1)+(SCw>1);
     SCp=SCp.*(SCp>=0).*(SCp<=1)+(SCp>1);
  i=j>=ndt; 
  
  iL=SCw_old(1:na)~=SCw(1:na);
    
 % i=1;
end;
%ndt