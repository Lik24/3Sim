function [X,A2X,dV,Won,L,gMp]=Rasp(Data,fl)
if fl==1
elseif fl==2
elseif fl==3
    X=Data.G;
    A2X=Data.A2G;
    dV=Data.dVg;
    %p=DData.pd;
    Won=Data.Won;
    L=Data.Lg;
    gMp=Data.gMp_g;    
else
    X=Data.D;
    A2X=Data.A2D;
    dV=Data.dVd;
    %p=DData.pd;
    Won=Data.Won;
    L=Data.Ld;
    gMp=Data.gMp_d;
end