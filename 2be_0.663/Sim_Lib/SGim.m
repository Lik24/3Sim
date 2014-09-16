function [Clp,Cwp,Cws,A,Bwo,Mp]=SGim(dV,SCw,dt,Mp,zc,Bwo,P,Patm,P0,va,vc,vg,vd)

Bw=Bwo(:,1:2);
Bo=Bwo(:,3:4);


[dVa,Sw,A_mp,BwA,BoA]=Narez(dV,SCw,Mp,Bw,Bo,va);
[dVc,Cw,C_mp,BwC,BoC]=Narez(dV,SCw,Mp,Bw,Bo,vc);
[dVg,Gw,G_mp,BwG,BoG]=Narez(dV,SCw,Mp,Bw,Bo,vg);
[dVd,Dw,D_mp,BwD,BoD]=Narez(dV,SCw,Mp,Bw,Bo,vd);

[CwpA,CopA,BwA,BoA,CwsA,AA,MpA]=ZaSgim(dVa,Sw,dt,A_mp,zc([1,2,4]),P(va),Patm,P0(va),BwA,BoA);
[CwpC,CopC,BwC,BoC,CwsC,AC,MpC]=ZaSgim(dVc,Cw,dt,C_mp,zc([1,2,4]),P(vc),Patm,P0(vc),BwC,BoC);
[CwpG,CopG,BwG,BoG,CwsG,AG,MpG]=ZaSgim(dVg,Gw,dt,G_mp,zc([1,2,4]),P(vg),Patm,P0(vg),BwG,BoG);
[CwpD,CopD,BwD,BoD,CwsD,AD,MpD]=ZaSgim(dVd,Dw,dt,D_mp,zc([1,2,6]),P(vd),Patm,P0(vd),BwD,BoD);

Bwo(:,1:2)=[BwA;BwC;BwG;BwD];
Bwo(:,3:4)=[BoA;BoC;BoG;BoD];

Cwp=[CwpA;CwpC;CwpG;CwpD];
Cop=[CopA;CopC;CopG;CopD];


Cws=[CwsA;CwsC;CwsG;CwsD];
A=[AA,AC,AG,AD];
Clp=A.*Cwp+Cop;

Mp=[MpA;MpC;MpG;MpD];

end

function [Cwp,Cop,Bw,Bo,Cws,A,Mp]=ZaSgim(dV,Sw,dt,mp,zc,P,Patm,P0,Bw,Bo)
    
    Bw2=Bw(:,1)./(1+zc(1)*(P-Patm));
    Bo2=Bo(:,1)./(1+zc(2)*(P-Patm));
    
    mp2=mp(:,1).*(1+zc(3)*(P-P0));
    
    Cwp=dV.*Sw/dt.*(mp(:,2)*zc(3)./Bw2+mp2*zc(1)./Bw(:,2));
    Cop=dV.*(1-Sw)/dt.*(mp(:,2)*zc(3)./Bo2+mp2*zc(2)./Bo(:,2));
    
    Cws=dV.*mp2./Bw2;
    A=Bw2./Bo2;
    
    Bw(:,2)=Bw2;
    Bo(:,2)=Bo2;
    Mp(:,1)=mp(:,1);
    Mp(:,2)=mp2;
end

function [dV,Sw,Mp,Bw,Bo]=Narez(dv,sw,mp,bw,bo,vec)

dV=dv(vec);
Sw=sw(vec);
Mp=mp(vec,:);
Bw=bw(vec,:);
Bo=bo(vec,:);
end