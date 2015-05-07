function [QQwoBND]=QIterGY2(dPt,QQwoBND,WBND,BXYZ)

QQwoBND.Qmw(BXYZ.mxy) = QQwoBND.Qmw(BXYZ.mxy) - WBND.b1gm(BXYZ.mxy,3).*dPt(BXYZ.mxy);
QQwoBND.Qmw(BXYZ.mz) = QQwoBND.Qmw(BXYZ.mz) - WBND.b1gm(BXYZ.mz,4).*dPt(BXYZ.mz);

QQwoBND.Qdw(BXYZ.dxy) = QQwoBND.Qdw(BXYZ.dxy) - WBND.b1gd(BXYZ.dxy,3).*dPt(BXYZ.dxy);
QQwoBND.Qdw(BXYZ.dz) = QQwoBND.Qdw(BXYZ.dz) - WBND.b1gd(BXYZ.dz,4).*dPt(BXYZ.dz);



