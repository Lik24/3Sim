function [G,A2G,dVg,p,WonG]=Gorizont(XY,SS,gt,WXY)

dh=0.0001;
[ntr,Nl]=size(gt);
Z=ones(size(XY,1),1);
Kg=1e-5;
np=size(XY,1);
[A]=MR_Prop(XY,1);
AS=A;
K=Kg*Z;      
WfG=zeros(size(WXY,1),ntr);
Wf=[];
R=[];
WNG=[];
Won=WfG;

for l=1:Nl
    GS=SS{l};
    HH=sparse(size(AS,1),size(AS,1));
    C2=sparse(size(AS,1),1);
    gt1=[];
    for i=1:ntr
        Gt=gt{i,l};
        H1=dh*ones(size(XY,1),1);
        v1=zeros(size(XY,1),1);
        v1(Gt)=1;
        nt=find(v1==0);
        
        [A]=MR_Prop(XY,1);
        [L,B,S,H]=Geome3(A,XY,Z,H1);
        % nt=randi(np,ceil(np*ch),1);
        if size(Gt)>3
            a=convhull(XY(Gt,1),XY(Gt,2));
            cxy=XY(Gt(a),:);
            Won(:,i)=inpolygon(WXY(:,1),WXY(:,2),cxy(:,1),cxy(:,2));
        else
            cxy(:,1)=0;
            cxy(:,2)=0;
            Won(:,i)=0;
        end;


        r=find(Won(:,i));
        WfG(r,i)=KWell(K,H1,sum(S,2),r);
        
        H(nt,:)=0;
        H(:,nt)=0;

        HH=H+HH;
        gt1=[gt1;Gt];
        C2=sum(GS,2);
 
        
       % figure(45),plot(cxy(:,1),cxy(:,2),WXY(:,1),WXY(:,2),'*')
       % hold on
    end;
    
        WONGM=sum(Won,2)~=0;
        r=find(WONGM);
        WFG=sum(WfG,2);
    
        Wf=[Wf;WFG(WONGM~=0)];
        R=[R;r+(l-1)*np];
        WNG=[WNG;r];
     
        v1=zeros(size(XY,1),1);
        v1(gt1)=1;
        nt=find(v1==0);
        C2(nt)=[];
            
        A(nt,:)=[];  A(:,nt)=[];
        L(nt,:)=[];  L(:,nt)=[];
        HH(nt,:)=[];  HH(:,nt)=[];
        B(nt,:)=[];  B(:,nt)=[];
        H=HH;
        
        v=1:np;
        v(nt)=[];
        
        n=size(A,1);
        [r,c]=find(A>0);
        
        G=H(r+(c-1)*n).*B(r+(c-1)*n)./L(r+(c-1)*n)*Kg;
        G=sparse(r,c,G,n,n);
        G=G-sparse(1:n,1:n,sum(G,2),n,n);
        
        bn=size(G,1);
        vDF=G(1:bn+1:end);
        c1=find(vDF==0);
        G(c1,:)=[];
        G(:,c1)=[];

        A2G=sparse(v,1:n,C2,np,n);
        A2G(:,c1)=[];
        dVg=sum(H.*L.*B,2);
        dVB(l)={dVg};
        CB(l)={G};
        A2CB(l)={A2G};
        
end;


dVg=[];
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
   dVg(size(dVg,2)+1:size(dVg,2)+ni)=dVB{i};
end;

G=cell2mat(cb);
A2G=cell2mat(a2cb);

    p=symrcm(G);
    G=G(p,p);
    A2G=A2G(:,p);
    dVg=dVg(p)';
    
for i=1:size(R,1)
  R(i)=find(R(i)==p);
end;

WonG(:,2)=Wf;
WonG(:,1)=R;
WonG(:,3)=WNG;
