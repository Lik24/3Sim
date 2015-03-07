function [Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,dt1,V0,Pw1,MJ,A1,Kxy]=Adap_GY(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,CR_GRUP,Pw_f)

dJ=0.1;
wxy=DATA.XY(DATA.Won,:);
for j=1:size(DATA.Won,1)
 R=((wxy(j,1)-wxy(:,1)).^2+(wxy(j,2)-wxy(:,2)).^2).^0.5;
 [r,I]=sort(R);
 blz(:,j)=I(2:5);
end;

[Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,dt1,V0]=SimT_MKT(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,1,CR_GRUP);
Pw1=Pw;
J1=-(Pw-Pw_f).*WData.Uf;
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
 
[A1,MJ,Kxy]=All_Well(mJ,PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,CR_GRUP,dJ,Pw_f,blz);
GYData.GY_Kxy=A1;

[Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,dt1,V0]=SimT_MKT(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,1,CR_GRUP);
end


function [AKX,MJ,Kxy]=All_Well(mJ,PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,CR_GRUP,dJ,Pw_f,blz)
  
 fl=max(abs(mJ))<dJ;
 k=0;
 mj=mJ;
 XY=DATA.XY;
 BndXY=DATA.BndXY;
 Kxy0=GYData.GY_Kxy(BndXY==1);
 
 while fl==0 && k<20
   k=k+1;
   
   for j=1:size(Pw_f,1)
    xy=DATA.XY(DATA.Won(j),:);
    R=((xy(1)-XY(BndXY==1,1)).^2+(xy(2)-XY(BndXY==1,2)).^2).^0.5;
    [~,I]=sort(R);
    
    Kxy(:,j)=GYData.GY_Kxy(BndXY==1);
    Kxy(I(1:5),j)=Kxy(I(1:5),j).*(1.15*(mJ(j)<-dJ)+0.85*(mJ(j)>dJ)+1*(abs(mJ(j))<=dJ));
    ZeroKxy=zeros(size(Kxy,1),1);
    ZeroKxy(I(1:5),1)=Kxy(I(1:5),j);
    Kxy(:,j)=ZeroKxy(:,1);
   end;
   
   for j=1:size(Kxy,1)
     v1=Kxy(j,:);
     KXY(j,1)=mean(Kxy(j,v1~=0));
     if isnan(KXY(j,1))==1
         KXY(j,1)=Kxy0(j);
     end;
   end;
%    plot(KXY)
%    hold on
   GYData.GY_Kxy(BndXY==1)=KXY;
      
   [Pi,Sw,Ti,MCp,p,Q,Pw,PpW,SwC,NDT,dQ,dSS,dt1,V0]=SimT_MKT(PR,C,A2C,G,A2G,dVc,dVg,DATA,WData,GYData,1,CR_GRUP);  
   J1=-(Pw-Pw_f).*WData.Uf;
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
 AKX=GYData.GY_Kxy;
 MJ=mj;
end