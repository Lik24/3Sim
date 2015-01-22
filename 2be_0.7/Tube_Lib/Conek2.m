function [C,A2C,dVc,p,WonV,L,CR_GRUP]=Conek2(XY,NLT,Nl,CrDATA,Won,r0,ka)

np=size(XY,1);
Z=ones(np,1);
HH=CrDATA.H;
dC=CrDATA.dC;
KC=CrDATA.KC;
DH=CrDATA.DH;
alp_C=CrDATA.alp_C;

WonV=zeros(0,3);

k=0;
nc1=0;
sntl=0;
[AA]=MR_Prop(XY,1);

dV_L=cell(Nl,1);
C_L=cell(Nl,1);
L_L=cell(Nl,1);
A2C_L=cell(Nl,1);

for l=1:Nl
    Nt=NLT{l};
    Kc=KC{l};
    Dh=DH{l};
    
    [L,~,~,H]=Geome3(AA,XY,Z,HH(:,l));
    mnt=size(Nt,2);
    dVB=cell(mnt,1);
    CB=cell(mnt,1);
    LB=cell(mnt,1);
    A2CB=cell(mnt,1);
    Cr_grup=cell(mnt,1);
    
    for i=1:mnt
        
        nt=Nt{i};   % Список ячеек с трещинами
        kc=Kc(i);   % Проницаемость трещины
        dh=Dh(i);   % Раскрытость трещины
        A=AA;
        A(:)=0;
         for ia=1:size(nt,2)
           A(nt(2,ia),nt(1,ia))=1;
         end;
         A=A+A';
%        figure(235),spy(A) 
       unt=unique(nt);

        A1=A(unt,:);   A2=A1(:,unt);
        L1=L(unt,:);   L2=L1(:,unt);
        H1=H(unt,:);   H2=H1(:,unt);
        
%        figure(234),spy(A2)
%        jhjh
%         [r,c]=find(A2==1);
%         for ir=1:size(r,1)
%            if sum(r(ir)==nt)==0
%             A2(r(ir),c(ir))=0;
%             L2(r(ir),c(ir))=0;
%             H2(r(ir),c(ir))=0;
%            end;
%         end;
%             unt'
%             Won'     
        for j=1:size(Won,1)
            ty=find(Won(j,1)==unt);
            if isempty(ty)==0
                k=k+1;
                WonV(k,1)=ty(1)+sntl;
                WonV(k,2)=HH(Won(j,1),l)*dh*kc*8.64/r0;%dC*100;
                WonV(k,3)=Won(j,3);
            end;
        end;
%         WonV
        n=size(A2,1);
        [r,c]=find(A2==1);
        C=H2(r+(c-1)*n)*dh./L2(r+(c-1)*n)*kc*8.34;
        C=sparse(r,c,C,n,n);
        
        a2c=sum(H2.*L2,2)*alp_C;
        A2C=sparse(unt+(l-1)*np,1:n,a2c,np*Nl,n);
        dVc=sum(H2.*L2.*dh,2);
        dVB(i)={dVc};
        CB(i)={C};
        LB(i)={L};
        A2CB(i)={A2C};
        Cr_grup(i)={[i*ones(size(C,1),1),l*ones(size(C,1),1)]};
        sntl=sntl+size(unt,1);
    end;

    nc1=nc1+sntl;
    C_L(l)=Mat_Constr(CB);
    L_L(l)=Mat_Constr(LB);
    A2C_L(l)={cell2mat(A2CB')};
    dV_L(l)={cell2mat(dVB)};
    CR_grup(l)={cell2mat(Cr_grup)};
end;

    C_cell=Mat_Constr(C_L);
    L_cell=Mat_Constr(L_L);
    A2C=cell2mat(A2C_L');
    dVc=cell2mat(dV_L);
    CR_GRUP=cell2mat(CR_grup');
    C=cell2mat(C_cell);
    L=cell2mat(L_cell);
   % A2C=cell2mat(A2C_cell);
   % dVc=cell2mat(dV_cell);


    p=symrcm(C);
    C=C(p,p);
    L=L(p,p);
    A2C=A2C(:,p);
    A2C=A2C(ka==1,:);
    dVc=dVc(p);
    CR_GRUP=CR_GRUP(p,:);
    for i=1:size(WonV,1)
        WonV(i,1)=find(WonV(i,1)==p);
    end;
end

function C=Mat_Constr(c)
n=size(c,1);
cb=cell(n,n);
for i=1:n
    A=c{i};
    ni=size(A,1);
    for j=1:n
        mj=size(c{j},2);
        cb(i,j)={sparse(ni,mj)};
    end;
    cb(i,i)=c(i);
end;

C=cell2mat(cb);
C={C};
end