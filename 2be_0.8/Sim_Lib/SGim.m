function [SGM]=SGim(dV,SCw,SCo,Mp,zc,Rs,rs,Bwog,dPt,P,Patm,va,vc,vg,vd,vb,dt)

Bw = Bwog(:,1:2);
Bo = Bwog(:,3:4);
Bg = Bwog(:,5:6);

[dVa,Sw,So,A_mp,RsA,BwA,BoA,BgA]=Narez(dV,SCw,SCo,Mp,Rs,Bw,Bo,Bg,va);
[dVc,Cw,Co,C_mp,RsC,BwC,BoC,BgC]=Narez(dV,SCw,SCo,Mp,Rs,Bw,Bo,Bg,vc);
[dVg,Gw,Go,G_mp,RsG,BwG,BoG,BgG]=Narez(dV,SCw,SCo,Mp,Rs,Bw,Bo,Bg,vg);
[dVd,Dw,Do,D_mp,RsD,BwD,BoD,BgD]=Narez(dV,SCw,SCo,Mp,Rs,Bw,Bo,Bg,vd);
[dVb,Bw,Bo,B_mp,RsB,BwB,BoB,BgB]=Narez(dV,SCw,SCo,Mp,Rs,Bw,Bo,Bg,vb);

[CwpA,CopA,CgpA,BwA,BoA,BgA,RsA,CgswA,CgsoA,CwsA,CosA,MpA]=ZaSgim(dVa,Sw,So,A_mp,zc([1,2,3,4]),RsA,rs,dPt(va),P(va),Patm,BwA,BoA,BgA,dt);
[CwpC,CopC,CgpC,BwC,BoC,BgC,RsC,CgswC,CgsoC,CwsC,CosC,MpC]=ZaSgim(dVc,Cw,Co,C_mp,zc([1,2,3,4]),RsC,rs,dPt(vc'),P(vc),Patm,BwC,BoC,BgC,dt);
[CwpG,CopG,CgpG,BwG,BoG,BgG,RsG,CgswG,CgsoG,CwsG,CosG,MpG]=ZaSgim(dVg,Gw,Go,G_mp,zc([1,2,3,4]),RsG,rs,dPt(vg'),P(vg),Patm,BwG,BoG,BgG,dt);
[CwpD,CopD,CgpD,BwD,BoD,BgD,RsD,CgswD,CgsoD,CwsD,CosD,MpD]=ZaSgim(dVd,Dw,Do,D_mp,zc([1,2,3,6]),RsD,rs,dPt(vd'),P(vd),Patm,BwD,BoD,BgD,dt);
[CwpB,CopB,CgpB,BwB,BoB,BgB,RsB,CgswB,CgsoB,CwsB,CosB,MpB]=ZaSgim(dVb,Bw,Bo,B_mp,zc([1,2,3,4]),RsB,rs,dPt(vb'),P(vb),Patm,BwB,BoB,BgB,dt);

SGM.Bwog(:,1:2)=[BwA;BwC;BwG;BwD;BwB];
SGM.Bwog(:,3:4)=[BoA;BoC;BoG;BoD;BoB];
SGM.Bwog(:,5:6)=[BgA;BgC;BgG;BgD;BgB];

SGM.Rs=[RsA;RsC;RsG;RsD;RsB];

SGM.Cwp=[CwpA;CwpC;CwpG;CwpD;CwpB];
SGM.Cop=[CopA;CopC;CopG;CopD;CopB];
SGM.Cgp=[CgpA;CgpC;CgpG;CgpD;CgpB];

SGM.Cgsw=[CgswA;CgswC;CgswG;CgswD;CgswB];
SGM.Cgso=[CgsoA;CgsoC;CgsoG;CgsoD;CgsoB];

SGM.Cws=[CwsA;CwsC;CwsG;CwsD;CwsB];
SGM.Cos=[CosA;CosC;CosG;CosD;CosB];

SGM.Clp = SGM.Cgp - SGM.Cgsw.*SGM.Cwp - SGM.Cgso.*SGM.Cop;

SGM.Mp=[MpA;MpC;MpG;MpD;MpB];
end

function [Cwp,Cop,Cgp,Bw,Bo,Bg,Rs,Cgsw,Cgso,Cws,Cos,Mp]=ZaSgim(dV,Sw,So,mp,zc,Rs,rs,dPt,P,Patm,Bw,Bo,Bg,dt)
    
    Bw2 = Bw(:,1).*(1 - zc(1)./(1 + zc(1)*(P - Patm)).*dPt);
    Bo2 = Bo(:,1).*(1 - zc(2)./(1 + zc(2)*(P - Patm)).*dPt);
    Bg2 = Bg(:,1).*(1 - zc(3)./(1 + zc(3)*(P - Patm)).*dPt);
    Rs2 = Rs + rs*dPt;
    mp2 = mp + zc(4)*dPt;
    
    Cwp = dV/dt.*Sw.*(mp*zc(4)./Bw2 + mp2*zc(1)./Bw(:,1));
    Cop = dV/dt.*So.*(mp*zc(4)./Bo2 + mp2*zc(2)./Bo(:,1));
    Cgp = dV/dt.*(mp*zc(4).*((1-Sw-So)./Bg2 + Rs2.*So./Bo2) + mp2.*(zc(3)./Bg(:,1) + So.*(rs./Bo2 + zc(2)*Rs2./Bo(:,1))));
    
    Cgsw = - Bw2./Bg2;
    Cgso = - Bo2./Bg2 + Rs2;
    
    Cws = dV/dt.*mp2./Bw2;
    Cos = dV/dt.*mp2./Bo2;
      
    Bw(:,2) = Bw2;
    Bo(:,2) = Bo2;
    Bg(:,2) = Bg2;
    Mp(:,1) = mp;
    Mp(:,2) = mp2;
    Rs(:,1) = Rs;
    Rs(:,2) = Rs2;
end

function [dV,Sw,So,Mp,Rs,Bw,Bo,Bg]=Narez(dv,sw,so,mp,rs,bw,bo,bg,vec)

dV=dv(vec);
Sw=sw(vec);
So=so(vec);
Mp=mp(vec);
Bw=bw(vec,:);
Bo=bo(vec,:);
Bg=bg(vec,:);
Rs=rs(vec);
end