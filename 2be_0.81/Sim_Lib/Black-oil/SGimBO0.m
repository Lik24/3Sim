function [CMP]=SGimBO0(Sw,So,zc,rs,P,Pb,CMP,V,dt,dV)
  
  zcb = 0.001;
    
  Rs(V.vp,1) = rs*Pb(V.vp);
  Rs(V.vg,1) = rs*P(V.vg);
  
  Bw = CMP.B0(1)./(1 + zc(1)*(P - Pb));
  Bo = (CMP.B0(2) + zcb*Rs)./(1 + zc(2)*(P - Pb));
  Bg = CMP.B0(3)./(1 + zc(3)*(P - Pb));
  Mp = CMP.Mp0.*(1 + zc(4)*(P - CMP.P0));
    
  dBoP = - CMP.B0(2)*zc(2)./(1 + zc(2)*(P - Pb));
  dBoPb = zcb*rs./(1 + zc(2)*(P - Pb)) - dBoP;  
  
  Coso = Mp./Bo;
  Cosw = - Coso;
  Cwsw = Mp./Bw;
   
  Cgsw(V.vg,1) = - Mp(V.vg)./Bg(V.vg);
  Cgso(V.vg,1) = - Bo(V.vg)./Bg(V.vg) + Rs(V.vg);
  Cgsw(V.vp,1) = - Mp(V.vp)./Bo(V.vp).*Rs(V.vp);
    
  Copb(V.vp,1) = -Mp(V.vp).*So(V.vp).*dBoPb(V.vp)./Bo(V.vp)./Bo(V.vp);
  Cgpb(V.vp,1) = Mp(V.vp).*So(V.vp).*(Bo(V.vp)*rs - Rs(V.vp).*dBoPb(V.vp))./Bo(V.vp)./Bo(V.vp);
  Db(V.vp,1) = Cosw(V.vp).*Cgpb(V.vp) - Cgsw(V.vp).*Copb(V.vp);
  
  Copb(V.vp) = Copb(V.vp)./Db(V.vp);
  Cgpb(V.vp) = Cgpb(V.vp)./Db(V.vp);
  
  CMP.Cg = ones(size(Sw,1),1);
  CMP.Cor(V.vp,1) = Copb(V.vp).*Cwsw(V.vp);
  CMP.Cor(V.vg,1) = CMP.Cg(V.vg);
  CMP.Co(V.vp,1) = -Cgpb(V.vp).*Cwsw(V.vp);
  CMP.Co(V.vg,1) = -Cgso(V.vg);
  CMP.Cw(V.vp,1) = CMP.Cg(V.vp);
  CMP.Cw(V.vg,1) = -Cgsw(V.vg)./Cwsw(V.vg);
  CMP.Cg(V.vp) = zeros(size(V.vp,1),1);
  
  CMP.Bw = repmat(Bw,1,2);
  CMP.Bo = repmat(Bo,1,2);
  CMP.Bg = repmat(Bg,1,2);
  CMP.Mp = repmat(Mp,1,2);
  CMP.Rs = repmat(Rs,1,2);
