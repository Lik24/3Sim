function [C,A2C,dVc,p,WonV]=Conek(XY,Nt,Nl,CrDATA,Won,dh,Kc,r0)

np=size(XY,1);
Z=ones(np,1);
H1=CrDATA.H;
dC=CrDATA.dC;
WonV(1,1:3)=1;
WonV(1,:)=[];

k=0;
nc1=0;

for i=1:Nl
    nt1=Nt{i}-np*(i-1);
    nt=1:np;
        
    nt(nt1)=[];
%     XY(nt1,:)'
    [A]=MR_Prop(XY,1);
    [L,B,S,H]=Geome3(A,XY,Z,H1(:,i));
    
    for j=1:size(Won,1)
        ty=find(Won(j)==nt1);
        if isempty(ty)==0
            k=k+1;
            WonV(k,1)=ty(1)+nc1;
            WonV(k,2)=H1(Won(j),i)*dh*Kc*8.64/r0;%dC*100;
            WonV(k,3)=Won(j);
        end;
    end;
    nc1=nc1+size(nt1,2);
    
    A2=A(nt1,:);
    A3=A2(:,nt1);

    L2=L(nt1,:);
    L3=L2(:,nt1);
    H2=H(nt1,:);
    H3=H2(:,nt1);
  %  full(A3)
    
    A(nt,:)=[];  A(:,nt)=[];
     
    L(nt,:)=[];  L(:,nt)=[];
    H(nt,:)=[];  H(:,nt)=[];
     
    A=A3;
    L=L3;
    H=H3;
    n=size(A,1);
    [r,c]=find(A==1);
% [r,c]'
    C=H(r+(c-1)*n)*dh./L(r+(c-1)*n)*Kc*8.34;
    C=sparse(r,c,C,n,n);
    C=C-sparse(1:n,1:n,sum(C,2),n,n);

    bn=size(C,1);
    vDF=C(1:bn+1:end);
    c1=find(vDF==0);
% c1
%     jhj
    C(c1,:)=[];
    C(:,c1)=[];
    
    C2=sum(H.*L,2);
    A2C=sparse(nt1,1:n,C2,np,n);
    A2C(:,c1)=[];
%     [r,c]=find(A2C);
%     [r,c]

    dVc=sum(H.*L.*dh,2);
    dVc(c1)=[];
    dVB(i)={dVc};
    CB(i)={C};
    %size(C)
    A2CB(i)={A2C};
end;
dVc=[];
for i=1:Nl
  [ni,mi]=size(CB{i});
  [nk,mk]=size(A2CB{i});
   for j=1:Nl
    [nj,mj]=size(CB{j});
    cb(i,j)={sparse(ni,mj)};
    [nj,mj]=size(A2CB{j});
    a2cb(i,j)={sparse(nk,mj)};
   end;
   cb(i,i)=CB(i);
   a2cb(i,i)=A2CB(i);

   dVc(size(dVc,2)+1:size(dVc,2)+ni)=dVB{i};
end;

C=cell2mat(cb);
A2C=cell2mat(a2cb);

    p=symrcm(C);
    C=C(p,p);
    A2C=A2C(:,p);
    dVc=dVc(p)';

%     XY(Nt{1},:)
%     sp=Nt{1};
%     XY(sp(WonV(:,1)),:)
%     fgklk

    for i=1:size(WonV,1)
        WonV(i,1)=find(WonV(i,1)==p);
    end;
%   WonV
