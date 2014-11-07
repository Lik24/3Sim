function [Doly,DATA,GYData]=Load_adp_prm(DATA,GYData,XY_new)

SD=load('ADP_ReZ2');
Doly(1:33)=SD.Doly;
gKX=SD.gKX;
Doly(34:36)=mean(SD.Doly./gKX(1:33)).*DATA.gKX(34:36);

GY_Kxy=SD.GY_Kxy(:,6);
GY_Pxy=SD.GY_Pxy;
P=SD.Pi0;
Sw=SD.Sw0;

P0=razmaz(P,SD.XY_old,XY_new);
Sw0=razmaz(Sw,SD.XY_old,XY_new);
gKX0=razmaz(gKX,SD.XY_old,XY_new);


DATA.gP=P0;
DATA.gSw=Sw0;
DATA.gKX=gKX0;

GY_Kxy=razmaz(GY_Kxy,SD.XY_old,XY_new);
GY_Pxy=razmaz(GY_Pxy,SD.XY_old,XY_new);

GYData.GY_Kxy=GY_Kxy;
GYData.GY_Pxy=GY_Pxy;
end

function B=razmaz(A,WXY,XY)
 F=scatteredInterpolant(WXY(:,1),WXY(:,2),A,'nearest','nearest');
 B=F(XY(:,1),XY(:,2));
end