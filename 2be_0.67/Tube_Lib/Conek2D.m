function [DData,CR_GRUP,Mp]=Conek2D(DATA,NLT,Nl,CrDATA,WData)

XY=DATA.XY;
Won=DATA.Won;
ka=DATA.ka;
BND=DATA.BND;
Hi=DATA.gH;
Z=DATA.gZ;
np=size(XY,1);
ddol=0.5;

KD=CrDATA.KD;
DH=CrDATA.DH;
Xk=CrDATA.Xk;

Nw=size(WData.Doly,1);
uj=repmat(1:Nw,1,Nl);

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
[L,B,S,H]=Geome3_1(AA,XY,Z,Hi);
ka(sum(AA)==-1)=0;
[r,c]=find(AA(1:size(XY,1),1:size(XY,1))==1);
K=cell2mat(KD);
K=repmat(K,size(XY,1),1);

K1=DATA.gKX*Xk;

[r1,c1]=find(AA==1);
KcKl=K1(r1)+K1(c1);
Ke=2*K1(r1).*K1(c1)./KcKl;
KK=sparse(r1,c1,Ke);

Wf=KWell(K1,ddol*Hi,S,L,B,Won,r,c,WData.Doly,WData.SDoly,WData.r0,XY,Nw,Nl);

im=zeros(size(ka));
im(Won)=uj;
im=im(ka==1);
uj=im(im~=0);

im=zeros(size(ka));
im(Won)=1;
im=im(ka==1);
Won=find(im(:));

for l=1:Nl
    Nt=NLT{l};
    Kd=KD{l};
    Dh=DH{l};
    
    mnt=size(Nt,2);
    dVB=cell(mnt,1);
    CB=cell(mnt,1);
    LB=cell(mnt,1);
    A2DB=cell(mnt,1);
    Cr_grup=cell(mnt,1);
    
    for i=1:mnt
        
        nt=Nt{i};   % Список ячеек с трещинами
        kd=Kd(i);   % Проницаемость трещины
        A=AA;
       
        unt=unique(nt);

        A1=A(unt,:);    A2=A1(:,unt);
        L1=L(unt,:);    L2=L1(:,unt);
        H1=H(unt,:);    H2=H1(:,unt);
        B1=B(unt,:);    B2=B1(:,unt);
        S1=S(unt,:);    S2=S1(:,unt);
        KK1=KK(unt,:);  KK2=KK1(:,unt);
    
        for j=1:size(Won,1)
            ty=find(Won(j)==unt);
            if isempty(ty)==0
                k=k+1;
                WonV(k,1)=ty(1)+sntl;
                WonV(k,2)=Wf(j);%dC*100;
                WonV(k,3)=uj(j);
            end;
        end;
%         WonV
        n=size(A2,1);
        [r,c]=find(A2==1);
        
        Lm=L2(r+(c-1)*n);
        Bm=B2(r+(c-1)*n);
        H1m=H2(r+(c-1)*n);
        K1m=KK2(r+(c-1)*n);
       % D=kd.*Bm.*H1m./Lm;
        km=K1m;
        D=km.*Bm.*H1m./Lm;
        %D=H2(r+(c-1)*n)*dh./L2(r+(c-1)*n)*kd*8.34;
        D=sparse(r,c,D,n,n);
        
        a2d=sum(H2.*S2,2)*1e-4;
        A2D=sparse(unt,1:n,a2d,np*Nl,n);
        dVd=sum(H2.*S2,2);
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
    
    if isempty(D)==0
        D=D(ka==1,ka==1);
        L=L(ka==1,ka==1);
        A2D=A2D(:,ka==1);
        dVd=dVd(ka==1);
    end
    
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
    
    Mp=DATA.gMp;
    Mp_d=ddol.*Mp(ka==1);
    Mp_d=Mp_d(p);
    Mp=ddol.*Mp;
       
    
DData.D=D;
DData.A2D=A2D;
DData.dVd=dVd;
DData.pd=p;
DData.Won=WonV;
DData.Ld=L;
DData.gMp_d=Mp_d;

DData.Kd(:,1)=DATA.gKX(ka==1)*Xk;
DData.Kd(:,2)=DATA.gKY(ka==1)*Xk;
DData.Kd(:,3)=DATA.gKZ(ka==1)*Xk;
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