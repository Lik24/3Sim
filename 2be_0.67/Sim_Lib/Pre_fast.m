function [TM,TC,TG,TD,TA2C,TA2G,TA2D,RC,Txyz_GY_A,Txyz_GY_D,dV1,dV2]=Pre_fast(A,C,G,D,A2C,A2G,A2D,C2G,Ke,L,B,S,H1,K,KeGY,BndXY,BndZ,nb,dV)
na=size(A,1); 
nc=size(C,1); 
ng=size(G,1); 
nd=size(D,1);

  [r,c]=find(A==1);  
  [r2,c2]=find(C>0);
  [r3,c3]=find(G>0);
  [r4,c4]=find(D>0);
  [r5,c5]=find(A2C);
  [r6,c6]=find(A2G);
  [r7,c7]=find(C2G);
  [r8,c8]=find(A2D);
     
  RC.rc(:,1)=[r;r2+na;c5+na;r5;r3+na+nc;c6+nc+na;r6];
  RC.rc(:,2)=[c;c2+na;r5;c5+na;c3+na+nc;r6;c6+nc+na];
  
Lm=L(r+(c-1)*na);
Bm=B(r+(c-1)*na);
H1m=H1(r+(c-1)*na);

TM=Ke.*Bm.*H1m./Lm;
RC.Arc=[r,c];
RC1=sparse(r,c,1);
U=triu(RC1);
[rh,ch]=find(U);
RC.Acr2=[rh,ch];

%TC=C(r2+(c2-1)*nc);
RC.Cr=r2;
RC.Cc=c2;

RC2=sparse(r2,c2,1);
U=triu(RC2);
[r2_h,c2_h]=find(U);
RC.Cr2=r2_h;
RC.Cc2=c2_h;
TC=C(r2_h+(c2_h-1)*nc);

%TG=G(r5+(c5-1)*ng);
RC.Gr=r3;
RC.Gc=c3;
RC5=sparse(r3,c3,1);
U=triu(RC5);
[r5_h,c5_h]=find(U);
RC.Gr2=r5_h;
RC.Gc2=c5_h;
TG=C(r5_h+(c5_h-1)*ng);

RC.Dr=r4;
RC.Dc=c4;
RC8=sparse(r4,c4,1);
U=triu(RC8);
[r8_h,c8_h]=find(U);
RC.Dr2=r8_h;
RC.Dc2=c8_h;
TD=D(r8_h+(c8_h-1)*nd);

TA2C=A2C(r5+(c5-1)*na).*K(r5);
RC.ACr=r5;
RC.ACc=c5;

TA2G=A2G(r6+(c6-1)*na).*K(r6);
RC.AGr=r6;
RC.AGc=c6;

TA2D=A2D(r8+(c8-1)*na).*K(r8);
RC.ADr=r8;
RC.ADc=c8;

TC2G=C2G(r6+(c6-1)*na).*K(r6);
RC.CGr=r7;
RC.CGc=c7;

%TA2G=A2G(r4+(c4-1)*na).*K(r4);
RC.ADr=r8;
RC.ADc=c8;

RC.na=na;
RC.nc=nc;
RC.ng=ng;
RC.nd=nd;
RC.nb=nb;

TM=full(TM);
TC=full(TC);
TG=full(TG);
TD=full(TD);

TA2C=full(TA2C);
TA2G=full(TA2G);
TA2D=full(TA2D);

vad=RC.ADr;
AGYXY=sparse(r,c,BndXY(c));
[rxy,cxy]=find(AGYXY);
BHL_XY=sparse(rxy,cxy,B(rxy+(cxy-1)*na).*H1(rxy+(cxy-1)*na)./L(rxy+(cxy-1)*na),na,na);


AGYZ=sparse(r,c,BndZ(c));
[rz,cz]=find(AGYZ);
BHLZ=sparse(rz,cz,B(rz+(cz-1)*na).*H1(rz+(cz-1)*na)./H1(rz+(cz-1)*na),na,na);
BlHXY=sum(BHL_XY,2);
BlHZ=sum(BHLZ,2);

Txyz_GY_A(:,1)=KeGY(:,1).*BlHXY;
Txyz_GY_A(:,2)=KeGY(:,2).*BlHZ;

Txyz_GY_D(:,1)=KeGY(vad,1).*BlHXY(vad);
Txyz_GY_D(:,2)=KeGY(vad,2).*BlHZ(vad);

dVa=dV(1:na);
dVd=dV(na+nc+ng+1:na+nc+ng+nd);

dV1(:,1)=dVa(RC.Acr2(:,2));
dV1(:,2)=dVa(RC.Acr2(:,1));
dV2(:,1)=dVd(RC.Dc2);
dV2(:,2)=dVd(RC.Dr2);