function [qwog]=QIterSeqBO(qwog,W6,Wo,Wg,dPt,Won,dPw,n)

dP = dPw - dPt(Won);
nw = size(dPw,1);
qwog(:,1) = qwog(:,1) + sparse(Won(:,1),ones(1,nw),W6.*dP,n,1);
qwog(:,2) = qwog(:,2) + sparse(Won(:,1),ones(1,nw),Wo.*dP,n,1);
qwog(:,3) = qwog(:,3) + sparse(Won(:,1),ones(1,nw),Wg.*dP,n,1);