function [QQ,qwo]=QIter(QQ,W6,Wo,Wg,dPi,Won,dPw,SGM)

dP=(dPw - dPi(Won));

qwo(:,1) = W6.*dP;
qwo(:,2) = Wo.*dP;
qg = Wg.*dP;

QQ = QQ + qg + (SGM.Rs(Won) - SGM.Cgso(Won)).*qwo(:,2) - SGM.Cgsw(Won).*qwo(:,1);