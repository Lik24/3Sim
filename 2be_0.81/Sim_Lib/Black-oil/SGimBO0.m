function [CMP]=SGimBO0(Sw,So,zc,rs,P,Pb,Patm,CMP,dV,dt,V)
  
  zcb = 0.001;
    
  Rs(V.vg,1) = rs*Pb(V.vg);
  Rs(V.vp,1) = rs*Pb(V.vp);
  
  Bw = CMP.Bo0(1)./(1 + zc(1)*(P - Pb));
  Bo = (CMP.Bo0(2) + zcb*Rs)./(1 + zc(2)*(P - Pb));
  Bg = CMP.Bo0(3)./(1 + zc(3)*(P - Pb));
  Mp = CMP.Mp0.*(1 + zc(4)*(P - CMP.P0));
  
  dBoP = - Bo*zc(2)./(1 + zc(2)*(P - Pb));
  dBoPb = zcb*rs./(1 + zc(2)*(P - Pb)) - dBoP;
    
  CMP.Bw = repmat(Bw,1,2);
  CMP.Bo = repmat(Bo,1,2);
  CMP.Bg = repmat(Bg,1,2);
  CMP.Mp = repmat(Mp,1,2);
  CMP.Rs = repmat(Rs,1,2);
