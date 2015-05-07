function [dItime,dItimeW] = NR_Time(dV,Sw,Sw0,CMP,dt)

dItimeW = dV./dt.*(CMP.Mp(:,2).*Sw./CMP.Bw(:,2) - CMP.Mp(:,1).*Sw0./CMP.Bw(:,1));
dItimeO = dV./dt.*(CMP.Mp(:,2).*(1-Sw)./CMP.Bo(:,2) - CMP.Mp(:,1).*(1-Sw0)./CMP.Bo(:,1));

dItime = CMP.Cw.*dItimeW + dItimeO;
