function [C,A2C]=Tube(Nl,XY,ch)

np=size(XY,1);
Z=ones(np,1);
H=ones(np,1);

DT=delaunayTriangulation(XY(:,1),XY(:,2));
[A]=MR_Prop(XY,1);
[L,B,S,H]=Geome3(A,XY,Z,H);

nt=randi(np,ceil(np*ch),1);
A(nt,:)=[];
A(:,nt)=[];
L(nt,:)=[];
L(:,nt)=[];
XY(nt,:)=[];
[r,c]=find(A);

[C,A2C]=Conek(A,L,H,nt,np,Nl);

nt=randi(size(A,1),1,1);
nt1=nt;
fl=zeros(size(XY,1),1);
fl(nt)=1;
fl1=fl;
k=0;
k1=0;
while k==0 && k1<100
 k1=k1+1;
 [r1,c1]=find(A(nt,:));
 c1(fl(c1)==1)=[];
 fl(c1)=1;
 nt=c1;
 fl2=fl;
 k=sum(fl1)==sum(fl2);
 fl1=fl2;
end;

[r2,c2]=find(A(fl1==1,fl1==1));
XY2=XY;
XY2(fl1~=1,:)=[];

At=A(fl1==1,fl1==1);
L(L~=0)=1./L(L~=0);
Lt=L(fl1==1,fl1==1);
n=sum(fl1);
[C,I1]=min(XY(fl1==1,1)+XY(fl1==1,2));
[C,I2]=max(XY(fl1==1,1)+XY(fl1==1,2));

Lt=Lt-sparse(1:n,1:n,sum(Lt),n,n)-sparse([I1,I2],[I1,I2],[2,2],n,n);
bt=zeros(sum(fl1),1);
bt(I1)=-10;
bt(I2)=-100;

p=bt'/Lt;
A2C=1;
VZL2be(XY,XY2,r2,c2,p,DT,I1,I2,nt1,fl1,r,c)