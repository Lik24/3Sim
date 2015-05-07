function [dItime,dItimeW] = NR_TimeF2(dV,Sw,Sw0,CMP,dt,va,vc,ndt)

dItimeW(va,1) = dV(va)./dt.*(CMP.Mp(va,2).*Sw(va)./CMP.Bw(va,2) - CMP.Mp(va,1).*Sw0(va)./CMP.Bw(va,1));
dItimeO(va,1) = dV(va)./dt.*(CMP.Mp(va,2).*(1-Sw(va))./CMP.Bo(va,2) - CMP.Mp(va,1).*(1-Sw0(va))./CMP.Bo(va,1));

dItimeW(vc,1) = ndt.*dV(vc)./dt.*(CMP.Mp(vc,2).*Sw(vc)./CMP.Bw(vc,2) - CMP.Mp(vc,1).*Sw0(vc)./CMP.Bw(vc,1));
dItimeO(vc,1) = ndt.*dV(vc)./dt.*(CMP.Mp(vc,2).*(1-Sw(vc))./CMP.Bo(vc,2) - CMP.Mp(vc,1).*(1-Sw0(vc))./CMP.Bo(vc,1));

dItime = CMP.Cw.*dItimeW + dItimeO;
