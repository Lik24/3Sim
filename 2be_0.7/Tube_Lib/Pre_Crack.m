function [CR_rc]=Pre_Crack(RC,na,nd,TM,TD,A2C,A2G,D2C,D2G,A2D,Wf,Won,WonM,WonD)
WonM=[Won,Wf,WonM];

 CR_rc(1,1)=conct2mat(na,RC.ACr,RC.AGr,RC.Arc2(:,2),RC.Arc2(:,1),TM,A2C,A2G,WonM);
 CR_rc(1,2)=conct2mat(nd,RC.DCc,RC.DGc,RC.Dc2,RC.Dr2,TD,D2C',D2G',WonD);
 
 v1=CR_rc(1,1).v;
 v2=CR_rc(1,2).v;
 a2d=A2D;
 a2d(v1==0,:)=[];
 a2d(:,v2==0)=[];
 [r,c,val]=find(a2d);
 CR_rc(1,3).r=r;
 CR_rc(1,3).c=c;
 CR_rc(1,3).a2d=val;
end

function [CR_rc,won]=conct2mat(n,inc,ing,r,c,T,A2C,A2G,WoM)
    rcm=unique([inc;ing]);
    v=zeros(n,1);
    v(rcm)=1;
    gh=ones(n,1);
    gh(rcm)=2;
    Img2=sparse(r,c,gh(c),n,n);
    Img2=Img2+Img2';

    [r1,c1]=find(Img2==2);
    de=[];
    for i=1:size(rcm,1)
        de=[de;find(rcm(i)==r1)];
    end;

    r_gy=r1;      r_gy(de,:)=[];     r_in=r1(de,:);
    c_gy=c1;      c_gy(de,:)=[];     c_in=c1(de,:);

    de=[];
    for i=1:size(rcm,1)
        de=[de;find(rcm(i)==c_in)];
    end;
    r_in=r_in(de);
    c_in=c_in(de);
    
    RC_IN=sparse(r_in,c_in,1);
    U=triu(RC_IN);
    [r_in_h,c_in_h]=find(U);

     T=sparse(r,c,T,n,n);
     T=T+T';
     
     T_gy=T(r_gy+n*(c_gy-1));
     T_in_h=T(r_in_h+n*(c_in_h-1));
     
     rc_gy=[r_gy,c_gy];
     rc_in_h=[r_in_h,c_in_h];
     
     A2C(v~=1,:)=[];
     [r1,c1]=find(A2C);

     A2G(v~=1,:)=[];
     [r2,c2]=find(A2G);
     
     qw=zeros(n,3);
     qw(WoM(:,1),1)=1;
     qw(WoM(:,1),2)=WoM(:,2);
     qw(WoM(:,1),3)=WoM(:,3);
     qw(v==0,:)=[];
     won(:,1)=find(qw(:,1));
     won(:,2)=qw(won(:,1),2);
     won(:,3)=qw(won(:,1),3);
    
     CR_rc.won=won;
     CR_rc.r1=r1;
     CR_rc.c1=c1;
     CR_rc.r2=r2;
     CR_rc.c2=c2;
     CR_rc.rc_gy=rc_gy;
     CR_rc.rc_in_h=rc_in_h;
     CR_rc.T_gy=T_gy;
     CR_rc.T_in_h=T_in_h;

     CR_rc.v=v;
end

function D2C=D_Conect(A2D,A2C)
 [r1,c1]=find(A2D);
 [r2,c2]=find(A2C);
 n1=size(A2C,2);
 n2=size(A2D,2);
 A1=sum(A2D,2);
 A2=sum(A2C,2);
 A=A1.*A2;
 r=find(A);
 R=[];
 C=[];
  for i=1:size(r,1)
      C(i)=c2(r(i)==r2);
      R(i)=c1(r(i)==r1);
  end
  D2C=sparse(C,R,1,n1,n2);
end