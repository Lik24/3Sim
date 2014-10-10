function CR=Narez_Index(v1,r1,na,r,c,T,Wf,Won,WonM,Wn,A2C,A2G)

v2=v1;
v2(r1)=1;
rcm=find(v2==1);
gh=zeros(na,1);
gh(rcm)=2;

Img2=sparse(r,c,gh(c),na,na);
[r1,c1]=find(Img2==2);
de=[];
for i=1:size(rcm,1)
    de=[de;find(rcm(i)==r1)];
end;
r_gy=r1;
c_gy=c1;
r_gy(de)=[];
c_gy(de)=[];

r_in=r1(de);
c_in=c1(de);

RC_IN=sparse(r_in,c_in,1);
U=triu(RC_IN);
[r_in_h,c_in_h]=find(U);

T=sparse(r,c,T);
T_gy=T(r_gy+na*(c_gy-1));
T_in_h=T(r_in_h+na*(c_in_h-1));

rc_gy=[r_gy,c_gy];
rc_in_h=[r_in_h,c_in_h];

A2C11=A2C;
A2C11(v2~=1,:)=[];
[r1,c1]=find(A2C11);

A2G11=A2G;
A2G11(v2~=1,:)=[];
[r2,c2]=find(A2G11);

qw=zeros(na,1);
qw(Won)=1;
qw(v1==0)=[];
won=find(qw);

qw=zeros(na,1);
qw(Won)=Wf;
qw(v2==0)=[];
wf=qw(won);

qw=zeros(na,1);
qw(Won)=WonM;
qw(v2==0)=[];
wn=qw(won);
hj=zeros(Wn,1);
hj(wn)=1;
wn1=find(hj==1);

CR.won=won;
CR.wf=wf;
CR.wn=wn;
CR.wn1=wn1;

CR.r1=r1;
CR.c1=c1;
CR.r2=r2;
CR.c2=c2;
CR.rc_gy=rc_gy;
CR.rc_in_h=rc_in_h;
CR.T_gy=T_gy;
CR.T_in_h=T_in_h;