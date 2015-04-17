function [Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,dt1,V0,Pw1,Pw2,MJ,AKX,A,MJ2]=Adap_Kswaz_All(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,CR_GRUP,Pw_f)

dJ=5;
[Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,dt1,V0]=SimT_MKT(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,1,CR_GRUP);
Pw1=Pw;
J1=(Pw-Pw_f).*WData.Uf;
for j=1:size(Pw,1)
 J2=J1(j,Pw_f(j,:)~=0);
 J3=J2(isnan(Pw_f(j,Pw_f(j,:)~=0))==0);
 mJ(j,1)=sum(J3);
end;
A=WData.Doly;

for fi=1:2
[A,MJ]=All_Well(mJ,PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,CR_GRUP,dJ,Pw_f);
WData.Doly=A;
[AKX,MJ2]=All_Well2(mJ,PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,CR_GRUP,dJ/2,Pw_f);
DATA.gKX(DATA.Won)=AKX;
end;
% WData.Doly=A;
% [Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,dt1,V0]=SimT_MKT(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,1,CR_GRUP);
 Pw2=Pw;
% J1=(Pw-Pw_f).*WData.Uf;
% for j=1:size(Pw,1)
%  J2=J1(j,Pw_f(j,:)~=0);
%  J3=J2(isnan(Pw_f(j,Pw_f(j,:)~=0))==0);
%  mJ(j,1)=sum(J3);
% end;
% 
% [AKX,MJ2]=All_Well2(mJ,PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,CR_GRUP,dJ/2,Pw_f);
% 
% 
% DATA.gKX(DATA.Won)=AKX;
[Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,dt1,V0]=SimT_MKT(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,1,CR_GRUP);
AKX=DATA.gKX;
end

function [A,MJ]=All_Well(mJ,PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,CR_GRUP,dJ,Pw_f)
 
 fl=max(abs(mJ))<dJ;
 k=0;
 mj=mJ;
 while fl==0 && k<2
   k=k+1;
   WData.Doly=WData.Doly.*(1.15*(mJ<-dJ)+0.85*(mJ>dJ)+1*(abs(mJ)<=dJ));
   WData.Doly((abs(mJ)<=dJ))=WData.Doly(abs(mJ)<=dJ).*(1.05*(mJ(abs(mJ)<=dJ)<-dJ/2)+0.95*(mJ(abs(mJ)<=dJ)>dJ/2)+1*(abs(mJ(abs(mJ)<=dJ))<=dJ/2));
   % WData.Doly(j)
   [Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,dt1,V0]=SimT_MKT(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,1,CR_GRUP);  
   J1=(Pw-Pw_f).*WData.Uf;
   for j=1:size(Pw,1)
       J2=J1(j,Pw_f(j,:)~=0);
       J3=J2(isnan(Pw_f(j,Pw_f(j,:)~=0))==0);
       mJ(j,1)=sum(J3);
   end;
   mj(:,k)=mJ;
 end
 A=WData.Doly;
 MJ={mj};
end

function [AKX,MJ]=All_Well2(mJ,PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,CR_GRUP,dJ,Pw_f)

 fl=max(abs(mJ))<dJ;
 k=0;
 mj=mJ;
 while fl==0 && k<2
   k=k+1;
   DATA.gKX(DATA.Won)=DATA.gKX(DATA.Won).*(1.15*(mJ<-dJ)+0.85*(mJ>dJ)+1*(abs(mJ)<=dJ));
  % WData.Doly((abs(mJ)<=dJ))=WData.Doly(abs(mJ)<=dJ).*(1.05*(mJ(abs(mJ)<=dJ)<-dJ/2)+0.95*(mJ(abs(mJ)<=dJ)>dJ/2)+1*(abs(mJ(abs(mJ)<=dJ))<=dJ/2));
   % WData.Doly(j)
   [Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,dt1,V0]=SimT_MKT(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,1,CR_GRUP);  
   J1=(Pw-Pw_f).*WData.Uf;
   for j=1:size(Pw,1)
       J2=J1(j,Pw_f(j,:)~=0);
       J3=J2(isnan(Pw_f(j,Pw_f(j,:)~=0))==0);
       mJ(j,1)=sum(J3);
   end;
   mj(:,k)=mJ;
 end
 AKX=DATA.gKX(DATA.Won);
 MJ={mj};
end