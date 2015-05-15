function [SGM,CMP]=SGimBO(dV,Sw,So,zc,rs,P,Pb,dPb,dPt,dt,V,CMP)
  
  zcb = 0.001;
  
  dBwP = - CMP.Bw(:,2)*zc(1)./(1 + zc(1)*(P - Pb));
  dBoP = - CMP.Bo(:,2)*zc(2)./(1 + zc(2)*(P - Pb));
  dBoPb = zcb*rs./(1 + zc(2)*(P - Pb)) - dBoP;
  dBgP = - CMP.Bg(:,2)*zc(3)./(1 + zc(3)*(P - Pb));

  Bw = CMP.Bw(:,2) + dBwP.*dPt;
  Bg = CMP.Bg(:,2) + dBgP.*dPt;
  Bo(V.vg,1) = CMP.Bo(V.vg,2) + dBoP(V.vg).*dPt(V.vg);
  Bo(V.vp,1) = CMP.Bo(V.vp,2) + dBoP(V.vp).*dPt(V.vp) + dBoPb(V.vp).*dPb(V.vp);
  
  Mp = CMP.Mp(:,2) + zc(4)*dPt;
 
  Rs(V.vg,1) = rs*P(V.vg);
  Rs(V.vp,1) = rs*Pb(V.vp);

% Rs(V.vg,1) = CMP.Rs(V.vg,2) + rs*dPt(V.vg);
% Rs(V.vp,1) = CMP.Rs(V.vp,2) + rs*dPb(V.vp);
   
  SGM.Cwp = dV/dt.*Sw.*(CMP.Mp0*zc(4)./Bw - Mp.*dBwP./Bw./Bw);
  SGM.Cop = dV/dt.*So.*(CMP.Mp0*zc(4)./Bo - Mp.*dBoP./Bo./Bo);
  SGM.Cgp(V.vg,1) = dV(V.vg)/dt.*(CMP.Mp0(V.vg)*zc(4).*((1-Sw(V.vg)-So(V.vg))./Bg(V.vg) + ...
  Rs(V.vg).*So(V.vg)./Bo(V.vg)) + Mp(V.vg).*(-(1-Sw(V.vg)-So(V.vg)).*dBgP(V.vg)./Bg(V.vg)./Bg(V.vg) + So(V.vg).*(-Rs(V.vg).*dBoP(V.vg)./Bo(V.vg)./Bo(V.vg) + rs./Bo(V.vg))));
  SGM.Cgp(V.vp,1) = dV(V.vp)/dt.*(CMP.Mp0(V.vp).*zc(4).*Rs(V.vp).*So(V.vp)./Bo(V.vp) - Mp(V.vp).*So(V.vp).*Rs(V.vp).*dBoP(V.vp)./Bo(V.vp)./Bo(V.vp));  
   
 
% SGM.Cwp = dV/dt.*Sw.*(CMP.Mp0*zc(4)./Bw + Mp.*zc(1)./CMP.B0(1));
 % SGM.Cop = dV/dt.*So.*(CMP.Mp0*zc(4)./Bo + Mp.*zc(2)./CMP.B0(2));
 % SGM.Cgp(V.vg,1) = dV(V.vg)/dt.*(CMP.Mp0(V.vg)*zc(4).*((1-Sw(V.vg)-So(V.vg))./Bg(V.vg) + ...
 % Rs(V.vg).*So(V.vg)./Bo(V.vg)) + Mp(V.vg).*((1-Sw(V.vg)-So(V.vg)).*zc(3)./CMP.B0(3) + So(V.vg).*(Rs(V.vg).*zc(2)./CMP.B0(2) + rs./Bo(V.vg))));
 % SGM.Cgp(V.vp,1) = dV(V.vp)/dt.*(CMP.Mp0(V.vp).*zc(4).*Rs(V.vp).*So(V.vp)./Bo(V.vp) + Mp(V.vp).*So(V.vp).*Rs(V.vp).*zc(2)./CMP.B0(2));  
  
  Coso = Mp./Bo;
  Cosw = - Coso;
  Cwsw = Mp./Bw;
    
  Cgsw(V.vg,1) = - Mp(V.vg)./Bg(V.vg);
  Cgso(V.vg,1) = - Bo(V.vg)./Bg(V.vg) + Rs(V.vg);
  Cgsw(V.vp,1) = - Mp(V.vp)./Bo(V.vp).*Rs(V.vp);
  Cgso(V.vp,1) = zeros(size(V.vp,1),1);

  Copb(V.vp,1) = -Mp(V.vp).*So(V.vp).*dBoPb(V.vp)./Bo(V.vp)./Bo(V.vp);
  Cgpb(V.vp,1) = Mp(V.vp).*So(V.vp).*(rs./Bo(V.vp) - Rs(V.vp).*dBoPb(V.vp)./Bo(V.vp)./Bo(V.vp));
 
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
  
  SGM.Clp(V.vg,1) = SGM.Cgp(V.vg) - Cgsw(V.vg)./Cwsw(V.vg).*SGM.Cwp(V.vg) - Cgso(V.vg).*SGM.Cop(V.vg);
  SGM.Clp(V.vp,1) = SGM.Cwp(V.vp) - Cwsw(V.vp).*Cgpb(V.vp).*SGM.Cop(V.vp) + Cwsw(V.vp).*Copb(V.vp).*SGM.Cgp(V.vp);
  
  SGM.Cwsw(V.vg,1) = dV(V.vg)./dt.*Cwsw(V.vg);
  SGM.Coso(V.vg,1) = dV(V.vg)./dt.*Coso(V.vg);
  
  SGM.Cosw(V.vp,1) = Cosw(V.vp)./Db(V.vp)./dV(V.vp).*dt;
  SGM.Cgsw(V.vp,1) = Cgsw(V.vp)./Db(V.vp)./dV(V.vp).*dt;
  
  SGM.Cgpb(V.vp,1) = Cgpb(V.vp)./dV(V.vp).*dt;
  SGM.Copb(V.vp,1) = Copb(V.vp)./dV(V.vp).*dt;
 
  CMP.Bw(:,2) = Bw;
  CMP.Bo(:,2) = Bo;
  CMP.Bg(:,2) = Bg;
  CMP.Mp(:,2) = Mp;
  CMP.Rs(:,2) = Rs;