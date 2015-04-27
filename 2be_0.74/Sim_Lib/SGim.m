function [Clp,Cwp,Cws,A,Bwo,Mp]=SGim(dV,Sw,Mp,zc,Bwo,P,Patm,P0,va,vc,vg,vd,vb,dt)

Bw=Bwo(:,1:2);
Bo=Bwo(:,3:4);

% [~,MSw,A_mp,BwA,BoA]=Narez(dV,Sw,Mp,Bw,Bo,va);
% [~,CSw,C_mp,BwC,BoC]=Narez(dV,Sw,Mp,Bw,Bo,vc);
% [~,GSw,G_mp,BwG,BoG]=Narez(dV,Sw,Mp,Bw,Bo,vg);
% [~,DSw,D_mp,BwD,BoD]=Narez(dV,Sw,Mp,Bw,Bo,vd);
% [~,BSw,B_mp,BwB,BoB]=Narez(dV,Sw,Mp,Bw,Bo,vb);

% [CwpA,CopA,BwA,BoA,CwsA,AA,MpA]=ZaSgim(MSw,A_mp,zc([1,2,4]),P(va),Patm,P0(va),BwA,BoA);
% [CwpC,CopC,BwC,BoC,CwsC,AC,MpC]=ZaSgim(CSw,C_mp,zc([1,2,4]),P(vc),Patm,P0(vc),BwC,BoC);
% [CwpG,CopG,BwG,BoG,CwsG,AG,MpG]=ZaSgim(GSw,G_mp,zc([1,2,4]),P(vg),Patm,P0(vg),BwG,BoG);
% [CwpD,CopD,BwD,BoD,CwsD,AD,MpD]=ZaSgim(DSw,D_mp,zc([1,2,6]),P(vd),Patm,P0(vd),BwD,BoD);
% [CwpB,CopB,BwB,BoB,CwsB,AB,MpB]=ZaSgim(BSw,B_mp,zc([1,2,4]),P(vb),Patm,P0(vb),BwB,BoB);

Zc([va,vc,vg,vd,vb],1)=zc(1);
Zc([va,vc,vg,vd,vb],2)=zc(2);
Zc([va,vc,vg,vd,vb],3)=zc(4);
Zc(vd,3)=zc(6);

[Cwp,Cop,Bw,Bo,Cws,A,Mp]=ZaSgim1(Sw,Mp,Zc,P([va,vc,vg,vd,vb]),Patm,P0,Bw,Bo);

% Bwo(:,1:2)=[BwA;BwC;BwG;BwD;BwB];
% Bwo(:,3:4)=[BoA;BoC;BoG;BoD;BoB];

% Cwp=dV.*[CwpA;CwpC;CwpG;CwpD;CwpB];
% Cop=dV.*[CopA;CopC;CopG;CopD;CopB];
% 
% Cws=dV.*[CwsA;CwsC;CwsG;CwsD;CwsB];

Bwo=[Bw,Bo];
Cwp=dV.*Cwp/dt;
Cop=dV.*Cop;
Cws=dV.*Cws;

%A=[AA;AC;AG;AD;AB];
Clp=A.*Cwp+Cop;
Clp=Clp/dt;

%Mp=[MpA;MpC;MpG;MpD;MpB];
end

function [Cwp,Cop,Bw,Bo,Cws,A,Mp]=ZaSgim(Sw,mp,zc,P,Patm,P0,Bw,Bo)
    
    Bw2=Bw(:,1)./(1+zc(1)*(P-Patm));
    Bo2=Bo(:,1)./(1+zc(2)*(P-Patm));
    
    mp2=mp(:,1).*(1+zc(3)*(P-P0));
    
    Cwp=Sw.*(mp(:,2)*zc(3)./Bw2+mp2*zc(1)./Bw(:,2));
    Cop=(1-Sw).*(mp(:,2)*zc(3)./Bo2+mp2*zc(2)./Bo(:,2));
    
    Cws=mp2./Bw2;
    A=Bw2./Bo2;
    
    Bw(:,2)=Bw2;
    Bo(:,2)=Bo2;
    Mp(:,1)=mp(:,1);
    Mp(:,2)=mp2;
end

function [Cwp,Cop,Bw,Bo,Cws,A,Mp]=ZaSgim1(Sw,mp,zc,P,Patm,P0,Bw,Bo)
    
    Bw2=Bw(:,1)./(1+zc(:,1).*(P-Patm));
    Bo2=Bo(:,1)./(1+zc(:,2).*(P-Patm));
    
    mp2=mp(:,1).*(1+zc(:,3).*(P-P0));
    
    Cwp=Sw.*(mp(:,2).*zc(:,3)./Bw2+mp2.*zc(:,1)./Bw(:,2));
    Cop=(1-Sw).*(mp(:,2).*zc(:,3)./Bo2+mp2.*zc(:,2)./Bo(:,2));
    
    Cws=mp2./Bw2;
    A=Bw2./Bo2;
    %A=Bo2./Bw2;
    
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