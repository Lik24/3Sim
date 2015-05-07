function [SGM,CMP]=SGimMD2(dV,Sw,zc,P,dPt,dt,CMP,v)
    
  dBoP = - CMP.Bo(v,2)*zc(2)./(1 + zc(2)*(P - CMP.P0(v))); 
  dBwP = - CMP.Bw(v,2)*zc(1)./(1 + zc(1)*(P - CMP.P0(v)));
 
  Bw = CMP.Bw(v,2) + dBwP.*dPt;
  Bo = CMP.Bo(v,2) + dBoP.*dPt; 
  Mp = CMP.Mp(v,2) + zc(4)*dPt;
   
  SGM.Cwp = dV/dt.*Sw.*(CMP.Mp0(v)*zc(4)./Bw - Mp.*dBwP./Bw./Bw);
  SGM.Cop = dV/dt.*(1-Sw).*(CMP.Mp0(v)*zc(4)./Bo - Mp.*dBoP./Bo./Bo); 
  
  CMP.Cw(v) = Bw./Bo;
  SGM.Clp = SGM.Cop + CMP.Cw(v).*SGM.Cwp;
  SGM.Cwsw = dV./dt.*Mp./Bw;
 
  CMP.Bw(v,2) = Bw;
  CMP.Bo(v,2) = Bo;
  CMP.Mp(v,2) = Mp;