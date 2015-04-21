dfunction [Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,dt1,V0,Pw1,Pw2,MJ,AKX,A]=Adap_Kswaz(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,CR_GRUP,Pw_f)

dJ=500;
[Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,dt1,V0]=SimT_MKT(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,1,CR_GRUP);
Pw1=Pw;
J1=(Pw-Pw_f).*WData.Uf;
for j=1:size(Pw,1)
 J2=J1(j,Pw_f(j,:)~=0);
 J3=J2(isnan(Pw_f(j,Pw_f(j,:)~=0))==0);
 mJ(j,1)=sum(J3);
end;
A=WData.Doly;
for j=33%1:size(Pw,1)
 [A(j),MJ(j)]=Single_Well(j,mJ(j),PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,CR_GRUP,dJ,Pw_f(j,:));
end;

WData.Doly=A;
[Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,dt1,V0]=SimT_MKT(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,1,CR_GRUP);
Pw2=Pw;
J1=(Pw-Pw_f).*WData.Uf;
for j=1:size(Pw,1)
 J2=J1(j,Pw_f(j,:)~=0);
 J3=J2(isnan(Pw_f(j,Pw_f(j,:)~=0))==0);
 mJ(j,1)=sum(J3);
end;
AKX=DATA.gKX(DATA.Won);
for j=33%1:size(Pw,1)
 [AKX(j),MJ(j)]=Single_Well2(j,mJ(j),PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,CR_GRUP,dJ/2,Pw_f(j,:));
end;

DATA.gKX(DATA.Won)=AKX;
[Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,dt1,V0]=SimT_MKT(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,1,CR_GRUP);

end

function [A,MJ]=Single_Well(j,mJ,PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,CR_GRUP,dJ,Pw_f)
 
 fl=abs(mJ)<dJ;
 k=0;
 mj=mJ;
 while fl==0 && k<10
   k=k+1;
   WData.Doly(j)=WData.Doly(j).*(1.15*(mJ<-dJ)+0.85*(mJ>dJ)+1*(abs(mJ)<=dJ));
   % WData.Doly(j)
   [Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,dt1,V0]=SimT_MKT(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,1,CR_GRUP);  
   J1=(Pw(j,:)-Pw_f).*WData.Uf(j,:);
   J2=J1(Pw_f~=0);
   J2=J2(isnan(J2)==0);
   mJ=sum(J2);
   mj(k)=mJ;
 end
 A=WData.Doly(j);
 MJ={mj};
end

function [AKX,MJ]=Single_Well2(j,mJ,PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,CR_GRUP,dJ,Pw_f)
 
 fl=abs(mJ)<dJ;
 k=0;
 mj=mJ;
 while fl==0 && k<1
   k=k+1;
   DATA.gKX(DATA.Won(j))=DATA.gKX(DATA.Won(j)).*(1.15*(mJ<-dJ)+0.85*(mJ>dJ)+1*(abs(mJ)<=dJ));
   % WData.Doly(j)
   [Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,dt1,V0]=SimT_MKT(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,1,CR_GRUP);  
   J1=(Pw(j,:)-Pw_f).*WData.Uf(j,:);
   J2=J1(Pw_f~=0);
   J2=J2(isnan(J2)==0);
   mJ=sum(J2);
   mj(k)=mJ;
 end
 AKX=DATA.gKX(DATA.Won(j));
 MJ={mj};
en