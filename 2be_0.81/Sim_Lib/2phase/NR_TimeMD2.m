function [dItime,dItimeW] = NR_TimeMD2(dV,Sw,Sw0,CMP,dt,v)

dItimeW = dV./dt.*(CMP.Mp(v,2).*Sw./CMP.Bw(v,2) - CMP.Mp(v,1).*Sw0./CMP.Bw(v,1));
dItimeO = dV./dt.*(CMP.Mp(v,2).*(1-Sw)./CMP.Bo(v,2) - CMP.Mp(v,1).*(1-Sw0)./CMP.Bo(v,1));

dItime = CMP.Cw(v).*dItimeW + dItimeO;
