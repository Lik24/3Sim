function [CR_rc]=Pre_Crack(RC,na,T,A2C,A2G,Wf,Won,WonM,Wn)

 v1=zeros(na,1);
 v1([RC.ACr;RC.AGr])=1;
 rcm=find(v1==1);
 
 r=RC.Arc(:,1);
 c=RC.Arc(:,2);

 gh=zeros(na,1);
 gh(rcm)=2;
 
 Img2=sparse(r,c,gh(c),na,na);
[r1,c1]=find(Img2==2);
% rcm
%[r1,c1]'
de=[];
 for i=1:size(rcm,1)
  de=[de;find(rcm(i)==r1)];
 end;
r_gy=r1;
c_gy=c1;
r_gy(de)=[];
c_gy(de)=[];
%[r_gy,c_gy]'

r_in=r1(de);
c_in=c1(de);

RC_IN=sparse(r_in,c_in,1);
U=triu(RC_IN);
[r_in_h,c_in_h]=find(U);

 T=sparse(r,c,T);
 T_gy=T(r_gy+na*(c_gy-1));
% T_in=T(r_in+na*(c_in-1));% Использование триуг матрицы вместо
 T_in_h=T(r_in_h+na*(c_in_h-1));
 
 rc_gy=[r_gy,c_gy];
% rc_in=[r_in,c_in]; % Использование триуг матрицы вместо
 rc_in_h=[r_in_h,c_in_h];
 
A2C11=A2C;
%A2C11=sparse(RC.ACr,RC.ACc,A2C11,na,nc);
A2C11(v1~=1,:)=[];
[r1,c1]=find(A2C11);

A2G11=A2G;
%A2G11=sparse(RC.AGr,RC.AGc,A2G11,na,ng);
A2G11(v1~=1,:)=[];
[r2,c2]=find(A2G11);

qw=zeros(na,1);
qw(Won)=1;
qw(v1==0)=[];
won=find(qw);

qw=zeros(na,1);
qw(Won)=Wf;
qw(v1==0)=[];
wf=qw(won);

qw=zeros(na,1);
qw(Won)=WonM;
qw(v1==0)=[];
wn=qw(won);
hj=zeros(Wn,1);
hj(wn)=1;
wn1=find(hj==1);

CR_rc.won=won;
CR_rc.wf=wf;
CR_rc.wn=wn;
CR_rc.wn1=wn1;

CR_rc.r1=r1;
CR_rc.c1=c1;
CR_rc.r2=r2;
CR_rc.c2=c2;
CR_rc.rc_gy=rc_gy;
%CR_rc.rc_in=rc_in;
CR_rc.rc_in_h=rc_in_h;
CR_rc.T_gy=T_gy;
%CR_rc.T_in=T_in;
CR_rc.T_in_h=T_in_h;