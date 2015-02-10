function [L,B,SV,H,HV]=Geome3_1(Ab,XY,Z,h)

X=XY(:,1);  Y=XY(:,2);
n=size(X,1);
Nl=size(Ab,1)/n;

A=Ab(1:n,1:n);
xc=sparse(n*n,1);
yc=xc; xl=xc; yl=xc; %B=xc;

[r,c]=find(A==1);
xl(r+(c-1)*n)=X(c);
yl(r+(c-1)*n)=Y(c);

xc(r+(c-1)*n)=X(r);
yc(r+(c-1)*n)=Y(r);

L=((xc-xl).^2+(yc-yl).^2).^0.5;
L=reshape(L,n,n);

A2=(A==1);
d=A2(r,:).*A2(c,:)~=0;
[r1,c1]=find(d);

xt=(X(r(r1))+X(c(r1))+X(c1))/3;
yt=(Y(r(r1))+Y(c(r1))+Y(c1))/3;

xf=(X(r(r1))+X(c(r1)))/2;
yf=(Y(r(r1))+Y(c(r1)))/2;

Con1=[xt,yt,c1;xf,yf,c(r1)];
%Con=contr(Con1);

a(:,1)=xt-xf;
a(:,2)=yt-yf;

xt(sum(a,2)==0)=[];
yt(sum(a,2)==0)=[];
xf(sum(a,2)==0)=[];
yf(sum(a,2)==0)=[];
r1(sum(a,2)==0)=[];
c1(sum(a,2)==0)=[];
a(sum(a,2)==0,:)=[];

b(:,1)=xl(r(r1)+(c(r1)-1)*n)-xc(r(r1)+(c(r1)-1)*n);
b(:,2)=yl(r(r1)+(c(r1)-1)*n)-yc(r(r1)+(c(r1)-1)*n);

% xe=X(c1);
% [xe(xt==xf),yt(xt==xf)]
% dfgh
LS=((xt-xf).^2+(yt-yf).^2).^0.5;

siny=abs(sum(a.*b,2)./(LS.*L(r(r1)+(c(r1)-1)*n)));
% siny(siny>1)-1
% dfgh
csn=(1-siny.^2).^0.5;
csn=real(csn);
dr=sparse(r1,c1,LS.*csn);
dr=sum(dr,2);
B=sparse(r,c,dr,n,n);

ss=polyarea([xc(r(r1)+(c(r1)-1)*n),xt,xf],[yc(r(r1)+(c(r1)-1)*n),yt,yf],2);
ds=sparse(r1,c1,ss);
ds=sum(ds,2);
S=sparse(r,c,ds,n,n);

Bc1=cell(Nl*Nl,1);
Bc1(:)={sparse(n,n)};
for i=1:Nl
    hc=h(r,i);
    hl=h(c,i);
    hm=(hc+hl)/2;
    dh=sparse(r,c,hm,n,n);
    Bc1(i+(i-1)*Nl)={dh};
end;
Bc1=reshape(Bc1,Nl,Nl);
H=cell2mat(Bc1);
HV=h(:);

L=A_constr(L,Nl,n);
B=A_constr(B,Nl,n);
S=A_constr(S,Nl,n);

v1=n*n*Nl+1:n*Nl+1:n*n*Nl*Nl;
v2=n+1:n*Nl+1:n*n*Nl*Nl-Nl*n*n;

% L(v1)=Z(n+1:end)-Z(1:end-n);
% L(v2)=Z(n+1:end)-Z(1:end-n);

SV=sum(S,2);

L=L+VRTL(Z(n+1:end)-Z(1:end-n),v1,v2,n*Nl);
B=B+VRTL(SV(1:(Nl-1)*n),v1,v2,n*Nl);
H=H+VRTL(ones(size(v1)),v1,v2,n*Nl);
end
function B=A_constr(A,Nl,ns)
Bc=cell(Nl*Nl,1);
Bc(:)={sparse(ns,ns)};
Bc(1:Nl+1:end)={A};
Bc=reshape(Bc,Nl,Nl);
B=cell2mat(Bc);
end
function B=VRTL(A,v1,v2,n)
B=sparse([v1,v2],ones(size([v1,v2])),[A,A],n*n,1);
B=reshape(B,n,n);
end
function [ Con1] = contr( B )

 q=zeros(size(B,1),1);
 
 for i=1:size(B,1) 
     for j=i+1:size(B,1) 
         if (B(i,1)==B(j,1)).*(B(i,2)==B(j,2)).*(B(i,3)==B(j,3)) 
             q(i)=1;
         end; 
     end;
 end;
 
 B1=B(:,1);
 B2=B(:,2); 
 B3=B(:,3); 
 B1(q==1)=[]; 
 B2(q==1)=[]; 
 B3(q==1)=[];
 clear B; 
 B=[B1,B2,B3];
 Con1=B;
end