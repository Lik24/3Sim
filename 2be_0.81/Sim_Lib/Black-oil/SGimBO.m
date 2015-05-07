function [SGM,CMP]=SGim(dV,Sw,So,zc,rs,P,Pb,dPb,dPt,dt,V,CMP)
  
  zcb = 0.001;
  
  dBoP = - CMP.Bo(:,2)*zc(2)./(1 + zc(2)*(P - Pb));
  dBoPb = zcb*rs./(1 + zc(2)*(P - Pb)) - dBoP;

  Bopb(V.vp,1) = CMP.Bo(V.vp,2);
  Bopb(V.vg,1) = CMP.Bo(V.vg,2) + dBoPb(V.vg).*dPb(V.vg);
 
  Bw = CMP.Bw(:,2).*(1 - zc(1)*dPt);
  Bg = CMP.Bg(:,2).*(1 - zc(3)*dPt);
  Bo = CMP.Bo(:,2) + dBoP.*dPt;
  
  Mp = CMP.Mp(:,2) + zc(4)*dPt;
  
  Rs(V.vp,1) = CMP.Rs(V.vp,2) + rs*dPt(V.vp);
  Rs(V.vg,1) = CMP.Rs(V.vg,2) + rs*dPb(V.vg);
   
  SGM.Cwp = dV/dt.*Sw.*(CMP.Mp0*zc(4)./Bw + Mp.*zc(1)./Bw);
  SGM.Cop = dV/dt.*So.*(CMP.Mp0*zc(4)./Bo + Mp.*zc(2)./Bo);
  SGM.Cgp(V.vp,1) = dV(V.vp)/dt.*(CMP.Mp0(V.vp)*zc(4).*((1-Sw(V.vp)-So(V.vp))./Bg(V.vp) + ...
  Rs(V.vp).*So(V.vp)./Bo(V.vp)) + Mp(V.vp).*((1-Sw(V.vp)-So(V.vp)).*zc(3)./Bg(V.vp) + So(V.vp).*Rs(V.vp).*zc(2)./Bo(V.vp)));
  SGM.Cgp(V.vg,1) = dV(V.vg)/dt.*(CMP.Mp0(V.vg).*zc(4).*Rs(V.vg).*So(V.vg)./Bo(V.vg) + Mp(V.vg).*So(V.vg).*Rs(V.vg).*zc(2)./Bo(V.vg));  
  
  Coso = Mp./Bo;
  Cosw = - Coso;
  Cwsw = Mp./Bw;
    
  Cgsw(V.vp,1) = - Mp(V.vp)./Bg(V.vp);
  Cgso(V.vp,1) = - Bo(V.vp)./Bg(V.vp) + Rs(V.vp);
  Cgsw(V.vg,1) = - Mp(V.vg)./Bo(V.vg).*Rs(V.vg);
  Cgso(V.vg,1) = zeros(size(V.vg,1),1);
 
  Copb(V.vg,1) = -Mp(V.vg).*So(V.vg).*(dBoPb(V.vg)./Bopb(V.vg)./CMP.Bo((V.vg),2));
  Cgpb(V.vg,1) = Mp(V.vg).*So(V.vg).*(CMP.Bo((V.vg),2)*rs - Rs(V.vg).*dBoPb(V.vg))./Bopb(V.vg)./CMP.Bo((V.vg),2);
  Db(V.vg,1) = Cosw(V.vg).*Cgpb(V.vg) - Cgsw(V.vg).*Copb(V.vg);
  
  Copb(V.vg) = Copb(V.vg)./Db(V.vg);
  Cgpb(V.vg) = Cgpb(V.vg)./Db(V.vg);
  
  SGM.Cg = ones(size(Sw,1),1);
  
  SGM.Cor(V.vg,1) = Copb(V.vg).*Cwsw(V.vg);
  SGM.Cor(V.vp,1) = SGM.Cg(V.vp);
  
  SGM.Co(V.vg,1) = -Cgpb(V.vg).*Cwsw(V.vg);
  SGM.Co(V.vp,1) = -Cgso(V.vp);

  SGM.Cw(V.vg,1) = SGM.Cg(V.vg);
  SGM.Cw(V.vp,1) = -Cgsw(V.vp)./Cwsw(V.vp);
  
  SGM.Cg(V.vg) = zeros(size(V.vg,1),1);
  
  SGM.Clp(V.vp,1) = SGM.Cgp(V.vp) - Cgsw(V.vp)./Cwsw(V.vp).*SGM.Cwp(V.vp) - Cgso(V.vp).*SGM.Cop(V.vp);
  SGM.Clp(V.vg,1) = SGM.Cwp(V.vg) - Cwsw(V.vg).*Cgpb(V.vg).*SGM.Cop(V.vg) + Cwsw(V.vg).*Copb(V.vg).*SGM.Cgp(V.vg);
  
  SGM.Cwsw(V.vp,1) = dV(V.vp)./dt.*Cwsw(V.vp);
  SGM.Coso(V.vp,1) = dV(V.vp)./dt.*Coso(V.vp);
  
  SGM.Cosw(V.vg,1) = Cosw(V.vg)./Db(V.vg)./dV(V.vg).*dt;
  SGM.Cgsw(V.vg,1) = Cgsw(V.vg)./Db(V.vg)./dV(V.vg).*dt;
  
  SGM.Cgpb(V.vg,1) = Cgpb(V.vg)./dV(V.vg).*dt;
  SGM.Copb(V.vg,1) = Copb(V.vg)./dV(V.vg).*dt;
 
  CMP.Bw(:,2) = Bw;
  CMP.Bo(:,2) = Bo + Bopb - CMP.Bo(:,2);
  CMP.Bg(:,2) = Bg;
  CMP.Mp(:,2) = Mp;
  CMP.Rs(:,2) = Rs;