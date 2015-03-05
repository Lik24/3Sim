function [SGM]=SGim0(SCw,SCo,zc,rs,P,Patm,P0,SGM,GEOM,VEC,dt)

Bw=SGM.Bwog(:,1:2);
Bo=SGM.Bwog(:,3:4);
Bg=SGM.Bwog(:,5:6);

[dVa,Sw,So,A_mp,BwA,BoA,BgA]=Narez(GEOM.dV,SCw,SCo,SGM.Mp,Bw,Bo,Bg,VEC.va);
[dVc,Cw,Co,C_mp,BwC,BoC,BgC]=Narez(GEOM.dV,SCw,SCo,SGM.Mp,Bw,Bo,Bg,VEC.vc);
[dVg,Gw,Go,G_mp,BwG,BoG,BgG]=Narez(GEOM.dV,SCw,SCo,SGM.Mp,Bw,Bo,Bg,VEC.vg);
[dVd,Dw,Do,D_mp,BwD,BoD,BgD]=Narez(GEOM.dV,SCw,SCo,SGM.Mp,Bw,Bo,Bg,VEC.vd);
[dVb,Bw,Bo,B_mp,BwB,BoB,BgB]=Narez(GEOM.dV,SCw,SCo,SGM.Mp,Bw,Bo,Bg,VEC.vb);

[CwpA,CopA,CgpA,BwA,BoA,BgA,RsA,CwsA,CosA,CgswA,CgsoA,MpA]=ZaSgim(dVa,Sw,So,A_mp,zc([1,2,3,4]),rs,P(VEC.va),Patm,P0(VEC.va),BwA,BoA,BgA,dt);
[CwpC,CopC,CgpC,BwC,BoC,BgC,RsC,CwsC,CosC,CgswC,CgsoC,MpC]=ZaSgim(dVc,Cw,Co,C_mp,zc([1,2,3,4]),rs,P(VEC.vc),Patm,P0(VEC.vc),BwC,BoC,BgC,dt);
[CwpG,CopG,CgpG,BwG,BoG,BgG,RsG,CwsG,CosG,CgswG,CgsoG,MpG]=ZaSgim(dVg,Gw,Go,G_mp,zc([1,2,3,4]),rs,P(VEC.vg),Patm,P0(VEC.vg),BwG,BoG,BgG,dt);
[CwpD,CopD,CgpD,BwD,BoD,BgD,RsD,CwsD,CosD,CgswD,CgsoD,MpD]=ZaSgim(dVd,Dw,Do,D_mp,zc([1,2,3,6]),rs,P(VEC.vd),Patm,P0(VEC.vd),BwD,BoD,BgD,dt);
[CwpB,CopB,CgpB,BwB,BoB,BgB,RsB,CwsB,CosB,CgswB,CgsoB,MpB]=ZaSgim(dVb,Bw,Bo,B_mp,zc([1,2,3,4]),rs,P(VEC.vb),Patm,P0(VEC.vb),BwB,BoB,BgB,dt);

SGM.Bwog(:,1:2)=[BwA;BwC;BwG;BwD;BwB];
SGM.Bwog(:,3:4)=[BoA;BoC;BoG;BoD;BoB];
SGM.Bwog(:,5:6)=[BgA;BgC;BgG;BgD;BgB];

SGM.Rs=[RsA;RsC;RsG;RsD;RsB];

SGM.Cwp=[CwpA;CwpC;CwpG;CwpD;CwpB];
SGM.Cop=[CopA;CopC;CopG;CopD;CopB];
SGM.Cgp=[CgpA;CgpC;CgpG;CgpD;CgpB];

SGM.Cws=[CwsA;CwsC;CwsG;CwsD;CwsB];
SGM.Cos=[CosA;CosC;CosG;CosD;CosB];

SGM.Cgsw=[CgswA;CgswC;CgswG;CgswD;CgswB];
SGM.Cgso=[CgsoA;CgsoC;CgsoG;CgsoD;CgsoB];

SGM.Clp = SGM.Cgp - SGM.Cgsw.*SGM.Cwp - SGM.Cgso.*SGM.Cop;

SGM.Mp=[MpA;MpC;MpG;MpD;MpB];
end

function [Cwp,Cop,Cgp,Bw,Bo,Bg,Rs,Cws,Cos,Cgsw,Cgso,Mp]=ZaSgim(dV,Sw,So,mp,zc,rs,P,Patm,P0,Bw,Bo,Bg,dt)
    
    Bw2 = Bw(:,1)./(1 + zc(1)*(P - Patm));
    Bo2 = Bo(:,1)./(1 + zc(2)*(P - Patm));
    Bg2 = Bg(:,1)./(1 + zc(3)*(P - Patm));
    Rs2 = rs*P;
    mp2 = mp(:,1).*(1 + zc(4)*(P - P0));
    
    Cwp = dV/dt.*Sw.*(mp(:,2)*zc(4)./Bw2 + mp2*zc(1)./Bw(:,2));
    Cop = dV/dt.*So.*(mp(:,2)*zc(4)./Bo2 + mp2*zc(2)./Bo(:,2));
    Cgp = dV/dt.*(mp(:,2)*zc(4).*((1-Sw-So)./Bg2 + Rs2.*So./Bo2) + mp2.*(zc(3)./Bg(:,2) + So.*(rs./Bo2 + zc(2)*Rs2./Bo(:,2))));
    
    Cos = dV/dt.*mp2./Bo2;
    Cws = dV/dt.*mp2./Bw2;
    
    Cgsw = - Bw2./Bg2;
    Cgso = - Bo2./Bg2 + Rs2;
      
    Bw(:,2)=Bw2;
    Bo(:,2)=Bo2;
    Bg(:,2)=Bg2;
    Mp(:,1)=mp(:,1);
    Mp(:,2)=mp2;
    Rs(:,2) = Rs2(:);
end

function [dV,Sw,So,Mp,Bw,Bo,Bg]=Narez(dv,sw,so,mp,bw,bo,bg,vec)

dV=dv(vec);
Sw=sw(vec);
So=so(vec);
Mp=mp(vec,:);
Bw=bw(vec,:);
Bo=bo(vec,:);
Bg=bg(vec,:);
end