function [QQwoBND]=QIterGY(dPt,QQwoBND,WBND,SGM,CMP,BXYZ)

QQwoBND.Qmw(BXYZ.mxy) = QQwoBND.Qmw(BXYZ.mxy) - WBND.b1gm(BXYZ.mxy,3).*dPt(BXYZ.mxy);
QQwoBND.Qmw(BXYZ.mz) = QQwoBND.Qmw(BXYZ.mz) - WBND.b1gm(BXYZ.mz,4).*dPt(BXYZ.mz);
QQwoBND.Qmo(BXYZ.mxy) = QQwoBND.Qmo(BXYZ.mxy) - WBND.b1gm(BXYZ.mxy,5).*dPt(BXYZ.mxy);
QQwoBND.Qmo(BXYZ.mz) = QQwoBND.Qmo(BXYZ.mz) - WBND.b1gm(BXYZ.mz,6).*dPt(BXYZ.mz);


QQwoBND.Qdw(BXYZ.dxy) = QQwoBND.Qdw(BXYZ.dxy) - WBND.b1gd(BXYZ.dxy,3).*dPt(BXYZ.dxy);
QQwoBND.Qdw(BXYZ.dz) = QQwoBND.Qdw(BXYZ.dz) - WBND.b1gd(BXYZ.dz,4).*dPt(BXYZ.dz);
QQwoBND.Qdo(BXYZ.dxy) = QQwoBND.Qdo(BXYZ.dxy) - WBND.b1gd(BXYZ.dxy,5).*dPt(BXYZ.dxy);
QQwoBND.Qdo(BXYZ.dz) = QQwoBND.Qdo(BXYZ.dz) - WBND.b1gd(BXYZ.dz,6).*dPt(BXYZ.dz);



