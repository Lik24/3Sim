function [vPa1,vPc1,vPg1,vPd1,dPa,dPc,dPg,dPd]=pre_potok_2(P,Pj,RC,rc_in_h,rc_in_hd,Na,Nd,na,nd)

%na=RC.na;
nc=RC.nc;
ng=RC.ng;

 r=RC.Arc(:,1);
 c=RC.Arc(:,2);
% 
%Pa=P(1:na);
% Pc=P(na+1:na+nc);
% Pg=P(na+nc+1:na+nc+ng);

dPa=P(rc_in_h(:,2),:)-P(rc_in_h(:,1),:);

vPa1(:,1)=dPa(:,1)>=0;
vPa1(:,2)=vPa1(:,1)==0;

vPa1(:,3)=dPa(:,2)>=0;
vPa1(:,4)=vPa1(:,3)==0;

dPd=P(Na+nc+ng+rc_in_hd(:,2),:)-P(Na+nc+ng+rc_in_hd(:,1),:);

vPd1(:,1)=dPd(:,1)>=0;
vPd1(:,2)=vPd1(:,1)==0;

vPd1(:,3)=dPd(:,2)>=0;
vPd1(:,4)=vPd1(:,3)==0;

if isempty(RC.Cr2)==1
    dPc=zeros(size(RC.Cr2,1),2);
    vPc1=zeros(size(RC.Cr2,1),4);
else
    dPc = Pj(na+RC.Cr2,:)-Pj(na+RC.Cc2,:);   
    vPc1(:,1)=dPc(:,1)>=0;
    vPc1(:,2)=vPc1(:,1)==0;
    vPc1(:,3)=dPc(:,2)>=0;
    vPc1(:,4)=vPc1(:,3)==0;
end

if isempty(RC.Gr2)==1
    dPg=zeros(size(RC.Gr2,1),2);
    vPg1=zeros(size(RC.Gr2,1),4);
else
    dPg = Pj(na+nc+RC.Gr2,:)-Pj(na+nc+RC.Gc2,:);
    vPg1(:,1)=dPg(:,1)>=0;
    vPg1(:,2)=vPg1(:,1)==0;
    vPg1(:,3)=dPg(:,2)>=0;
    vPg1(:,4)=vPg1(:,3)==0;
end