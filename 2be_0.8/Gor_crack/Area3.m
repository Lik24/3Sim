function [Z]=Area3(Ab,XY,cxy)

X=XY(:,1);  Y=XY(:,2);
n=size(X,1);

A=Ab(1:n,1:n);
xc=sparse(n*n,1);
yc=xc; %B=xc;

[r,c]=find(A==1);
xc(r+(c-1)*n)=X(r);
yc(r+(c-1)*n)=Y(r);


A2=(A==1);
d=A2(r,:).*A2(c,:)~=0;
[r1,c1]=find(d);

xt=(X(r(r1))+X(c(r1))+X(c1))/3;
yt=(Y(r(r1))+Y(c(r1))+Y(c1))/3;

xf=(X(r(r1))+X(c(r1)))/2;
yf=(Y(r(r1))+Y(c(r1)))/2;


ss=polyarea([xc(r(r1)+(c(r1)-1)*n),xt,xf],[yc(r(r1)+(c(r1)-1)*n),yt,yf],2);
ds=sparse(r1,c1,ss);
ds=sum(ds,2);
S=sparse(r,c,ds);

tx1=xc(r(r1)+(c(r1)-1)*n);
tx2=xt;
tx3=xf;

ty1=yc(r(r1)+(c(r1)-1)*n);
ty2=yt;
ty3=yf;

on1=inpolygon(tx1,ty1,cxy(:,1),cxy(:,2));
on2=inpolygon(tx2,ty2,cxy(:,1),cxy(:,2));
on3=inpolygon(tx3,ty3,cxy(:,1),cxy(:,2));

on=on1+on2+on3;
s1=(on==1);
s2=(on==2);
s3=(on==3);


ons1=sparse(r1,c1,s1);
ons2=sparse(r1,c1,s2);
ons3=sparse(r1,c1,s3);

ds1=sum(ons1,2);
ds2=sum(ons2,2);
ds3=sum(ons3,2);

mS1=sparse(r,c,ds1);
mS2=sparse(r,c,ds2);
mS3=sparse(r,c,ds3);

mS1=1/6*(mS1==1)+2/6*(mS1==2);
mS2=3/6*(mS2==1)+4/6*(mS2==2);
mS3=5/6*(mS3==1)+5/6*(mS3==2);

SS=S.*(mS1+mS2+mS3);
Z=sum(SS,2);