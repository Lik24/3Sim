function [Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,dt1,V0,Pw1,MJ,AKX,A]=Adap_Kswaz_All_2(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,CR_GRUP,Pw_f)

dJ=1;
wxy=DATA.XY(DATA.Won,:);
for j=1:size(DATA.Won,1)
 R=((wxy(j,1)-wxy(:,1)).^2+(wxy(j,2)-wxy(:,2)).^2).^0.5;
 [r,I]=sort(R);
 blz(:,j)=I(2:5);
end;

[Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,dt1,V0]=SimT_MKT(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,1,CR_GRUP);
Pw1=Pw;
J1=(Pw-Pw_f).*WData.Uf;
v1=[];
for j=1:size(Pw,1)
    J2=J1(j,Pw_f(j,:)~=0);
    J3=J2(isnan(Pw_f(j,Pw_f(j,:)~=0))==0);
    if isempty(J3)==0
        mJ(j,1)=sum(J3)/size(J3,2);
    else
       v1=[v1,j];
    end
end;
mJ(v1,1)=-sum(mJ(blz(:,v1)))/size(mJ(blz(:,v1)),2); 
 
A=WData.Doly;
MJ(:,1)=mJ;

for fi=1:10
[A,MJfi]=All_Well(MJ(:,end),PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,CR_GRUP,dJ,Pw_f,blz);
WData.Doly=A;
MJ=[MJ,MJfi];
[AKX,MJ2fi]=All_Well2(MJ(:,end),PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,CR_GRUP,dJ/2,Pw_f,blz);
DATA.gKX(DATA.Won)=AKX;
MJ=[MJ,MJ2fi];
end;

[Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,dt1,V0]=SimT_MKT(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,1,CR_GRUP);
AKX=DATA.gKX;
end

function [A,MJ]=All_Well(mJ,PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,CR_GRUP,dJ,Pw_f,blz)
 
 fl=max(abs(mJ))<dJ;
 k=0;
 mj=mJ;
 while fl==0 && k<2
   k=k+1;
   WData.Doly=WData.Doly.*(1.15*(mJ<-dJ)+0.85*(mJ>dJ)+1*(abs(mJ)<=dJ));
   % WData.Doly(j)
   [Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,dt1,V0]=SimT_MKT(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,1,CR_GRUP);  
   J1=(Pw-Pw_f).*WData.Uf;
   for j=1:size(Pw,1)
       J2=J1(j,Pw_f(j,:)~=0);
       J3=J2(isnan(Pw_f(j,Pw_f(j,:)~=0))==0);
       if isempty(J3)==0
        mJ(j,1)=sum(J3)/size(J3,2);
       else
         mJ(j,1)=-sum(mJ(blz(:,j)))/size(mJ(blz(:,j)),2); 
       end
   end;
   mj(:,k)=mJ;
 end
 A=WData.Doly;
 MJ=mj;
end

function [AKX,MJ]=All_Well2(mJ,PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,CR_GRUP,dJ,Pw_f,blz)

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
       if isempty(J3)==0
           mJ(j,1)=sum(J3)/size(J3,2);
       else
           mJ(j,1)=-sum(mJ(blz(:,j)))/size(mJ(blz(:,j)),2); 
       end
   end;
   mj(:,k)=mJ;
 end
 AKX=DATA.gKX(DATA.Won);
 MJ=mj;
end