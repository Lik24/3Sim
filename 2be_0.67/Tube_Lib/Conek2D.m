function [D,A2D,dVd,p,WonV,L,CR_GRUP]=Conek2D(DATA,NLT,Nl,CrDATA,r0)

XY=DATA.XY;
Won=DATA.Won;
ka=DATA.ka;
BND=DATA.BND;
H=DATA.gH;
Z=DATA.gZ;

np=size(XY,1);
HH=CrDATA.H;
dC=CrDATA.dC;
KD=CrDATA.KD;
DH=CrDATA.DH;


WonV(1,1:3)=1;
WonV(1,:)=[];

k=0;
nc1=0;
sntl=0;

dV_L=cell(Nl,1);
C_L=cell(Nl,1);
L_L=cell(Nl,1);
A2D_L=cell(Nl,1);

[AA]=MR_Prop_Bond(XY,Nl,BND);
[L,B,S,H1]=Geome3_1(AA,XY,Z,H);
ka(sum(AA)==-1)=0;
[r,c]=find(AA(1:size(XY,1),1:size(XY,1))==1);
Wf=KWell(KD,H,S,L,B,Won,r,c,WData.Doly,WData.r0,XY,Nw,Nl);

for l=1:Nl
    Nt=NLT{l};
    Kd=KD{l};
    Dh=DH{l};
    
    [L,~,~,H]=Geome3(AA,XY,Z,HH(:,l));
    mnt=size(Nt,2);
    dVB=cell(mnt,1);
    CB=cell(mnt,1);
    LB=cell(mnt,1);
    A2DB=cell(mnt,1);
    Cr_grup=cell(mnt,1);
    
    for i=1:mnt
        
        nt=Nt{i};   % Список ячеек с трещинами
        kd=Kd(i);   % Проницаемость трещины
        dh=Dh(i);   % Раскрытость трещины
        A=AA;
       
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
            ty=find(Won(j)==unt);
            if isempty(ty)==0
                k=k+1;
                WonV(k,1)=ty(1)+sntl;
                WonV(k,2)=Wf(Won(j));%dC*100;
                WonV(k,3)=Won(j);
            end;
        end;
%         WonV
        n=size(A2,1);
        [r,c]=find(A2==1);
        D=H2(r+(c-1)*n)*dh./L2(r+(c-1)*n)*kd*8.34;
        D=sparse(r,c,D,n,n);
        
        a2d=sum(H2.*L2,2);
        A2D=sparse(unt+(l-1)*np,1:n,a2d,np*Nl,n);
        dVd=sum(H2.*L2.*dh,2);
        dVB(i)={dVd};
        CB(i)={D};
        LB(i)={L};
        A2DB(i)={A2D};
        Cr_grup(i)={[i*ones(size(D,1),1),l*ones(size(D,1),1)]};
        sntl=sntl+size(unt,1);
    end;

    nc1=nc1+sntl;
    C_L(l)=Mat_Constr(CB);
    L_L(l)=Mat_Constr(LB);
    A2D_L(l)={cell2mat(A2DB')};
    dV_L(l)={cell2mat(dVB)};
    CR_grup(l)={cell2mat(Cr_grup)};
end;

    C_cell=Mat_Constr(C_L);
    L_cell=Mat_Constr(L_L);
    A2D=cell2mat(A2D_L');
    dVd=cell2mat(dV_L);
    CR_GRUP=cell2mat(CR_grup');
    D=cell2mat(C_cell);
    L=cell2mat(L_cell);
   % A2C=cell2mat(A2C_cell);
   % dVc=cell2mat(dV_cell);


    p=symrcm(D);
    D=D(p,p);
    L=L(p,p);
    A2D=A2D(:,p);
    A2D=A2D(ka==1,:);
    dVd=dVd(p);
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