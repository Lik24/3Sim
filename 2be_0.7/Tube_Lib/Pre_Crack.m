function [CR_rc]=Pre_Crack(RC,na,T,TD,A2C,A2G,Wf,Won,WonM,Wn)
WonM=[Won,Wf,WonM];

 [A,won]=conct2mat(na,RC.ACr,RC.AGr,RC.Arc(:,1),RC.Arc(:,2),T,A2C,A2G,WonM);
 
CR_rc.won=won;
CR_rc.wf=wf;
CR_rc.wn=wn;
CR_rc.wn1=wn1;

CR_rc.r1=r1;
CR_rc.c1=c1;
CR_rc.r2=r2;
CR_rc.c2=c2;
CR_rc.rc_gy=rc_gy;
CR_rc.rc_in=rc_in;
CR_rc.rc_in_h=rc_in_h;
CR_rc.T_gy=T_gy;
%CR_rc.T_in=T_in;
CR_rc.T_in_h=T_in_h;
CR_rc.TD_in_h=TD_in_h;
end

function [A,won]=conct2mat(n,inc,ing,r,c,T,A2C,A2G,WoM)
    rcm=unique([inc;ing]);
    v=zeros(n,1);
    v(rcm)=1;
    gh=zeros(n,1);
    gh(rcm)=2;
    Img2=sparse(r,c,gh(c),n,n);

    [r1,c1]=find(Img2==2);
    de=[];
    for i=1:size(rcm,1)
        de=[de;find(rcm(i)==r1)];
    end;

    r_gy=r1;        r_gy(de)=[];     r_in=r1(de);
    c_gy=c1;        c_gy(de)=[];     c_in=c1(de);
  
    RC_IN=sparse(r_in,c_in,1);
    U=triu(RC_IN);
    [r_in_h,c_in_h]=find(U);

     T=sparse(r,c,T);
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
    
end