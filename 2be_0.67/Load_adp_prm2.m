function [Doly,gKX,GY_Kxy]=Load_adp_prm2(Doly0,gy_Kxy,GY_Kxy)
SD=load('ADP_ReZ2');
Doly=SD.Doly;
gKX=SD.gKX;
GY_Kxy=SD.GY_Kxy(:,6);

