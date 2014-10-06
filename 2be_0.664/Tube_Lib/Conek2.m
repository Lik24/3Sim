function [C,A2C,dVc,p,WonV,L]=Conek2(XY,NLT,Nl,CrDATA,Won,r0)

np=size(XY,1);
Z=ones(np,1);
HH=CrDATA.H;
dC=CrDATA.dC;
KC=CrDATA.KC;
DH=CrDATA.DH;

WonV(1,1:3)=1;
WonV(1,:)=[];

k=0;
nc1=0;
sntl=0;
[A]=MR_Prop(XY,1);

dV_L=cell(Nl,1);
C_L=cell(Nl,1);
L_L=cell(Nl,1);
A2C_L=cell(Nl,1);

for l=1:Nl
    Nt=NLT{l};
    Kc=KC{l};
    Dh=DH{l};
    
    [L,~,~,H]=Geome3(A,XY,Z,HH(:,l));
    mnt=size(Nt,2);
    for i=1:mnt
        
        dVB=cell(mnt,1);
        CB=cell(mnt,1);
        LB=cell(mnt,1);
        A2CB=cell(mnt,1);

        nt=Nt{i};   % Список ячеек с трещинами
        kc=Kc(i);   % Проницаемость трещины
        dh=Dh(i);   % Раскрытость трещины
        
        A1=A(nt,:);   A2=A1(:,nt);
        L1=L(nt,:);   L2=L1(:,nt);
        H1=H(nt,:);   H2=H1(:,nt);
        
        for j=1:size(Won,1)
            ty=find(Won(j)==nt);
            if isempty(ty)==0
                k=k+1;
                WonV(k,1)=ty(1)+nc1;
                WonV(k,2)=HH(Won(j),l)*dh*kc*8.64/r0;%dC*100;
                WonV(k,3)=Won(j);
            end;
        end;
        
        n=size(A2,1);
        [r,c]=find(A2==1);
        C=H2(r+(c-1)*n)*dh./L2(r+(c-1)*n)*kc*8.34;
        C=sparse(r,c,C,n,n);
        
        a2c=sum(H2.*L2,2);
        A2C=sparse(nt,1:n,a2c,np,n);
        
        dVc=sum(H2.*L2.*dh,2);
        dVB(i)={dVc};
        CB(i)={C};
        LB(i)={L};
        A2CB(i)={A2C};
        
        sntl=sntl+size(nt,2);
    end;
    nc1=nc1+sntl;
    
    C_L(l)=Mat_Constr(CB);
    L_L(l)=Mat_Constr(LB);
    A2C_L(l)=Mat_Constr(A2CB);
    dV_L(l)=Mat_Constr(dVB);
end;

    C_cell=Mat_Constr(C_L);
    L_cell=Mat_Constr(L_L);
    A2C_cell=Mat_Constr(A2C_L);
    dV_cell=Mat_Constr(dV_L);
    C_cell
    sdsd
    C=cell2mat(C_cell);
    L=cell2mat(L_cell);
    A2C=cell2mat(A2C_cell);
    dVc=cell2mat(dV_cell);

    p=symrcm(C);
    C=C(p,p);
    L=L(p,p);
    A2C=A2C(:,p);
    dVc=dVc(p);

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

% for i=1:Nl
%     nt1=Nt{i}-np*(i-1);
%     nt=1:np;
% 
%     nt(nt1)=[];
% %     XY(nt1,:)'
%     
%     
%     
%     for j=1:size(Won,1)
%         ty=find(Won(j)==nt1);
%         if isempty(ty)==0
%             k=k+1;
%             WonV(k,1)=ty(1)+nc1;
%             WonV(k,2)=H1(Won(j),i)*dh*Kc*8.64/r0;%dC*100;
%             WonV(k,3)=Won(j);
%         end;
%     end;
%     nc1=nc1+size(nt1,2);
%     
%     A2=A(nt1,:);
%     A3=A2(:,nt1);
% 
%     L2=L(nt1,:);
%     L3=L2(:,nt1);
%     H2=H(nt1,:);
%     H3=H2(:,nt1);
%   %  full(A3)
%     
%     A(nt,:)=[];  A(:,nt)=[];
%      
%     L(nt,:)=[];  L(:,nt)=[];
%     H(nt,:)=[];  H(:,nt)=[];
%      
%     A=A3;
%     L=L3;
%     H=H3;
%     n=size(A,1);
%     [r,c]=find(A==1);
% % [r,c]'
%     C=H(r+(c-1)*n)*dh./L(r+(c-1)*n)*Kc*8.34;
%     C=sparse(r,c,C,n,n);
%     C=C-sparse(1:n,1:n,sum(C,2),n,n);
% 
%     bn=size(C,1);
%     vDF=C(1:bn+1:end);
%     c1=find(vDF==0);
% % c1
% %     jhj
%     C(c1,:)=[];
%     C(:,c1)=[];
%     
% 
%     kljgkj
%     L(c1,:)=[];
%     L(:,c1)=[];
%     
%     C2=sum(H.*L,2);
%     A2C=sparse(nt1,1:n,C2,np,n);
%     A2C(:,c1)=[];
% %     [r,c]=find(A2C);
% %     [r,c]
% 
%     dVc=sum(H.*L.*dh,2);
%     dVc(c1)=[];
%     dVB(i)={dVc};
%     CB(i)={C};
%     LB(i)={L};
%     %size(C)
%     A2CB(i)={A2C};
% end;
% dVc=[];
% for i=1:Nl
%   [ni,mi]=size(CB{i});
%   [nk,mk]=size(A2CB{i});
%    for j=1:Nl
%     [nj,mj]=size(CB{j});
%     cb(i,j)={sparse(ni,mj)};
%     lb(i,j)={sparse(ni,mj)};
%     [nj,mj]=size(A2CB{j});
%     a2cb(i,j)={sparse(nk,mj)};
%    end;
%    cb(i,i)=CB(i);
%    a2cb(i,i)=A2CB(i);
%    lb(i,i)=LB(i);
%    
%    dVc(size(dVc,2)+1:size(dVc,2)+ni)=dVB{i};
% end;
% 
% C=cell2mat(cb);
% L=cell2mat(lb);
% A2C=cell2mat(a2cb);
% 
%     p=symrcm(C);
%     C=C(p,p);
%     L=L(p,p);
%     A2C=A2C(:,p);
%     dVc=dVc(p)';
% 
% %     XY(Nt{1},:)
% %     sp=Nt{1};
% %     XY(sp(WonV(:,1)),:)
% %     fgklk
% 
%     for i=1:size(WonV,1)
%         WonV(i,1)=find(WonV(i,1)==p);
%     end;
% %   WonV
