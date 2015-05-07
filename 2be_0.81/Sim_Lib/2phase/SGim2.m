function [SGM,CMP]=SGim2(dV,Sw,zc,P,dPt,dt,CMP)
    
  dBoP = - CMP.Bo(:,2)*zc(2)./(1 + zc(2)*(P - CMP.P0)); 
  dBwP = - CMP.Bw(:,2)*zc(1)./(1 + zc(1)*(P - CMP.P0));
 
  Bw = CMP.Bw(:,2) + dBwP.*dPt;
  Bo = CMP.Bo(:,2) + dBoP.*dPt; 
  Mp = CMP.Mp(:,2) + zc(4)*dPt;
   
  SGM.Cwp = dV/dt.*Sw.*(CMP.Mp0*zc(4)./Bw - Mp.*dBwP./Bw./Bw);
  SGM.Cop = dV/dt.*(1-Sw).*(CMP.Mp0*zc(4)./Bo - Mp.*dBoP./Bo./Bo); 
  
  CMP.Cw = Bw./Bo;
  SGM.Clp = CMP.Cw.*SGM.Cwp + SGM.Cop;
  SGM.Cwsw = dV./dt.*Mp./Bw;
 
  CMP.Bw(:,2) = Bw;
  CMP.Bo(:,2) = Bo;
  CMP.Mp(:,2) = Mp;