function [CR_rc]=Pre_Crack(RC,na,nd,TM,TD,A2C,A2G,D2C,D2G,A2D,Wf,Won,WonM,WonD,dZ)
WonM=[Won,Wf,WonM];

 CR_rc(1,1)=conct2mat(na,RC.ACr,RC.AGr,RC.Arc(:,2),RC.Arc(:,1),TM,A2C,A2G,WonM,RC.Arc2(:,2),RC.Arc2(:,1));
 CR_rc(1,2)=conct2mat(nd,RC.DCr,RC.DGr,RC.Dc2,RC.Dr2,TD,D2C,D2G,WonD,RC.Dc,RC.Dr);
 
 v1=CR_rc(1,1).v;
 v2=CR_rc(1,2).v;
 a2d=A2D;
 a2d(v1==0,:)=[];
 a2d(:,v2==0)=[];
 [r,c,val]=find(a2d);
 CR_rc(1,3).r=r;
 CR_rc(1,3).c=c;
 CR_rc(1,3).a2d=val;
 
 dz=dZ(1,:);  dzw=dz{1};     dzo=dz{2};  dzg=dz{3};
 
 dzw=dzw(:,v1==1);   dzw=dzw(v1==1,:);
 dzo=dzo(:,v1==1);   dzo=dzo(v1==1,:);
 dzg=dzg(:,v1==1);   dzg=dzg(v1==1,:);
 
 dz(1)={dzw};  dz(2)={dzo};  dz(3)={dzg};  dZ(1,:)=dz;
   
 dz=dZ(4,:);  dzw=dz{1};     dzo=dz{2};  dzg=dz{3};
 
 dzw=dzw(:,v2==1);   dzw=dzw(v2==1,:);
 dzo=dzo(:,v2==1);   dzo=dzo(v2==1,:);
 dzg=dzg(:,v2==1);   dzg=dzg(v2==1,:);
 
 dz(1)={dzw};  dz(2)={dzo};  dz(3)={dzg};  dZ(4,:)=dz;
 CR_rc(1,4).dZ=dZ;
end

function [CR_rc,won]=conct2mat(n,inc,ing,c,r,T,A2C,A2G,WoM,ch,rh)
 if isempty(r)==0
    rcm=unique([inc;ing]);
    
    [C,ai,bi]=intersect(rcm,WoM(:,1));
    if isempty(C)==0
        coWi=WoM(bi,3);
        uCow=unique(coWi);
        inw=[];
        for i=uCow'
            inw=[inw;WoM(i==WoM(:,3),1)];
        end
        rcm=unique([rcm;inw]);
    end
    v=zeros(n,1);
    v(rcm)=1;

    Img2=sparse(r,c,v(c),n,n);
    %Img2=Img2+Img2';
    [r1,c1]=find(Img2==1);
    de=[];
    for i=1:size(rcm,1)
        de=[de;find(rcm(i)==r1)];
    end;
    
    [~,ci,de2]=intersect(rcm,r1,'rows','stable');
    [la,lb]=ismember(r1,rcm);
    
    
    r_gy=r1;     % r_gy(de3,:)=[];     r_in=r1(de3,:);
    c_gy=c1;     % c_gy(de3,:)=[];     c_in=c1(de3,:);

        r_gy=r1(la==0);
        c_gy=c1(la==0);
        
        r_in=r1(la==1);
        c_in=c1(la==1);     
        
%     de=[];
%     for i=1:size(rcm,1)
%         de=[de;find(rcm(i)==c_in)];
%     end;
%     [~,ci,de1]=intersect(rcm,c_in,'rows','stable');
%     r_in=r_in(de1);
%     c_in=c_in(de1);
[la,lb]=ismember(c_in,rcm);
        r_in=r_in(la==1);
        c_in=c_in(la==1);   
    
    RC_IN=sparse(r_in,c_in,1,n,n);
    U=tril(RC_IN);
    [r_in_h,c_in_h]=find(U);

     T=sparse(rh,ch,T,n,n);
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
 else
     won=zeros(0,3);
     CR_rc.won=won;
     CR_rc.r1=zeros(0,1);
     CR_rc.c1=zeros(0,1);
     CR_rc.r2=zeros(0,1);
     CR_rc.c2=zeros(0,1);
     CR_rc.rc_gy=zeros(0,2);
     CR_rc.rc_in_h=zeros(0,2);
     CR_rc.T_gy=zeros(0,1);
     CR_rc.T_in_h=zeros(0,1);
     CR_rc.v=zeros(0,1);
 end
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