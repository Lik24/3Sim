function [dItime,dItimeW,dItimeO,dItimeRs] = NR_Time(dV,Sw,So,Sw0,So0,SGM,CMP,dt)

dItimeW = dV./dt.*(CMP.Mp(:,2).*Sw./CMP.Bw(:,2) - CMP.Mp(:,1).*Sw0./CMP.Bw(:,1));
dItimeO = dV./dt.*(CMP.Mp(:,2).*So./CMP.Bo(:,2) - CMP.Mp(:,1).*So0./CMP.Bo(:,1));
dItimeRs = dV./dt.*(CMP.Mp(:,2).*CMP.Rs(:,2).*So./CMP.Bo(:,2) - CMP.Mp(:,1).*CMP.Rs(:,1).*So0./CMP.Bo(:,1));
dItimeG = dV./dt.*(CMP.Mp(:,2).*(1-Sw-So)./CMP.Bg(:,2) - CMP.Mp(:,1).*(1-Sw0-So0)./CMP.Bg(:,1));

dItime = SGM.Cg.*dItimeG + SGM.Cw.*dItimeW + SGM.Cor.*dItimeRs + SGM.Co.*dItimeO;
