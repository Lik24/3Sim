function [SGM,CMP]=SGimF2(dV,Sw,zc,P,dPt,dt,CMP,va,vc,ndt)
    
  dBoP = - CMP.Bo(:,2)*zc(2)./(1 + zc(2)*(P - CMP.P0)); 
  dBwP = - CMP.Bw(:,2)*zc(1)./(1 + zc(1)*(P - CMP.P0));
 
  Bw = CMP.Bw(:,2) + dBwP.*dPt;
  Bo = CMP.Bo(:,2) + dBoP.*dPt; 
  Mp = CMP.Mp(:,2) + zc(4)*dPt;
   
  SGM.Cwp(va,1) = dV(va)/dt.*Sw(va).*(CMP.Mp0(va)*zc(4)./Bw(va) - Mp(va).*dBwP(va)./Bw(va)./Bw(va));
  SGM.Cop(va,1) = dV(va)/dt.*(1-Sw(va)).*(CMP.Mp0(va)*zc(4)./Bo(va) - Mp(va).*dBoP(va)./Bo(va)./Bo(va)); 
  
  SGM.Cwp(vc,1) = dV(vc).*ndt./dt.*Sw(vc).*(CMP.Mp0(vc)*zc(4)./Bw(vc) - Mp(vc).*dBwP(vc)./Bw(vc)./Bw(vc));
  SGM.Cop(vc,1) = dV(vc).*ndt./dt.*(1-Sw(vc)).*(CMP.Mp0(vc)*zc(4)./Bo(vc) - Mp(vc).*dBoP(vc)./Bo(vc)./Bo(vc));
  
  CMP.Cw = Bw./Bo;
  SGM.Clp = CMP.Cw.*SGM.Cwp + SGM.Cop;
  SGM.Cwsw(va,1) = dV(va)./dt.*Mp(va)./Bw(va);
  SGM.Cwsw(vc,1) = dV(vc).*ndt./dt.*Mp(vc)./Bw(vc);
 
  CMP.Bw(:,2) = Bw;
  CMP.Bo(:,2) = Bo;
  CMP.Mp(:,2) = Mp;