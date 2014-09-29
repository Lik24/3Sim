function [NB,PXY]=Tube_perc(PR,CrDATA,XY,ch,WXY)

Nl=PR.Nl;
H=CrDATA.H;

np=size(XY,1);
Z=ones(np,1);

DT=delaunayTriangulation(XY(:,1),XY(:,2));
X_Y=XY;

for ij=1:Nl
    XY=X_Y;
    [A]=MR_Prop(XY,1);
    L=Geome3(A,XY,Z,H(:,ij));
 %   ij
 %   if ij==1
        nb=randi(np,ceil(np*ch),1);
        %sd=(XY(:,2)<160).*(XY(:,2)>080).*(XY(:,1)<960).*(XY(:,1)>040);
        %nb=find(sd(:)==0);
        NB(ij)={nb+np*(ij-1)};
%     else
%         sd=(XY(:,1)<160).*(XY(:,1)>080).*(XY(:,2)<960).*(XY(:,2)>040);
%         nb=find(sd(:)==0);
%         NB(ij)={nb};
%     end
%     
    A(nb,:)=[];
    A(:,nb)=[];
    L(nb,:)=[];
    L(:,nb)=[];
    XY(nb,:)=[];
    [r,c]=find(A);
    

    nt=randi(size(A,1),1,1);
    nt1=nt;
    fl=zeros(size(XY,1),1);
    fl(nt)=1;
    fl1=fl;
    k=0;
    k1=0;
    nt_nw=[];
    while k==0 && k1<100
        k1=k1+1;
        
        for i=1:size(nt,2)
            c1=find(A(nt(i),:));
            c1(fl(c1)==1)=[];
            fl(c1)=1;
            if isempty(nt_nw)==1
                nt_nw=c1;
            else
                nt_nw(size(nt_nw,2)+1:size(nt_nw,2)+size(c1,2))=c1;
            end;
        end;
        nt=nt_nw;
        nt_nw=[];
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

    PXY(:,ij)=VZL2be(XY,XY2,r2,c2,p,DT,I1,I2,nt1,fl1,r,c);
   
end;

for i=1:Nl
    v=1:np;
    n=NB{i}-np*(i-1);
    v(n)=[];
    NB(i)={v+np*(i-1)};
end;