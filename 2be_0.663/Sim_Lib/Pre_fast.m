function [TM,TC,TG,TA2C,TA2G,RC,T_GY]=Pre_fast(A,C,A2C,A2G,G,Ke,L,B,H1,K,KeGY,bz,rz,cz)
na=size(A,1); 
nc=size(C,1); 
ng=size(G,1); 

  [r,c]=find(A==1);  
  [r2,c2]=find(C>0);
  [r3,c3]=find(A2C);
  [r4,c4]=find(A2G);
  [r5,c5]=find(G>0);
     
  RC.rc(:,1)=[r;r2+na;c3+na;r3;r5+na+nc;c4+nc+na;r4];
  RC.rc(:,2)=[c;c2+na;r3;c3+na;c5+na+nc;r4;c4+nc+na];
  
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
RC.Gr=r5;
RC.Gc=c5;
RC5=sparse(r5,c5,1);
U=triu(RC5);
[r5_h,c5_h]=find(U);
RC.Gr2=r5_h;
RC.Gc2=c5_h;
TG=C(r5_h+(c5_h-1)*ng);

TA2C=A2C(r3+(c3-1)*na).*K(r3);
RC.ACr=r3;
RC.ACc=c3;

TA2G=A2G(r4+(c4-1)*na).*K(r4);
RC.AGr=r4;
RC.AGc=c4;

RC.na=na;
RC.nc=nc;
RC.ng=ng;

TM=full(TM);
TC=full(TC);
TG=full(TG);
TA2C=full(TA2C);
TA2G=full(TA2G);

BHLZ=sparse(rz,cz,B(rz+(cz-1)*na).*H1(rz+(cz-1)*na)./L(rz+(cz-1)*na),na,na);

T_GY=KeGY.*sum(BHLZ,2);