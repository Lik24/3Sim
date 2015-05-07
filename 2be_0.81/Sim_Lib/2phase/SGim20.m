function [CMP]=SGim20(Sw,So,zc,P,CMP,dV,dt)
     
  Bw = CMP.Bo0(1)./(1 + zc(1)*(P - CMP.P0));
  Bo = CMP.Bo0(2)./(1 + zc(2)*(P - CMP.P0));
  Mp = CMP.Mp0.*(1 + zc(4)*(P - CMP.P0));
    
  CMP.Bw = repmat(Bw,1,2);
  CMP.Bo = repmat(Bo,1,2);
  CMP.Mp = repmat(Mp,1,2);
  CMP.Cw = Bw./Bo;