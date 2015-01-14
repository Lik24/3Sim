function [vPa1,vPc1,vPg1,vPd1,dPa,dPc,dPg,dPd]=pre_potok_2(P,Pj,RC,rc_in_h,rc_in_hd,Na,Nd,na,nd,dZ)

%na=RC.na;
nc=RC.nc;
ng=RC.ng;

 r=RC.Arc(:,1);
 c=RC.Arc(:,2);
% 
%Pa=P(1:na);
% Pc=P(na+1:na+nc);
% Pg=P(na+nc+1:na+nc+ng);

dPa=P(rc_in_h(:,2))-P(rc_in_h(:,1));
DZ=dZ(1,:);
dZW=DZ{1};  dZO=DZ{2};  dZG=DZ{3};
vPa1(:,1)=dPa+dZW(rc_in_h(:,1)+(rc_in_h(:,2)-1)*na)>=0;
vPa1(:,2)=vPa1(:,1)==0;

vPa1(:,3)=dPa+dZO(rc_in_h(:,1)+(rc_in_h(:,2)-1)*na)>=0;
vPa1(:,4)=vPa1(:,3)==0;

dPd=P(Na+nc+ng+rc_in_hd(:,2))-P(Na+nc+ng+rc_in_hd(:,1));
DZ=dZ(4,:);
dZW=DZ{1};  dZO=DZ{2};  dZG=DZ{3};
vPd1(:,1)=dPd+dZW(rc_in_hd(:,1)+(rc_in_hd(:,2)-1)*nd)>=0;
vPd1(:,2)=vPd1(:,1)==0;

vPd1(:,3)=dPd+dZO(rc_in_hd(:,1)+(rc_in_hd(:,2)-1)*nd)>=0;
vPd1(:,4)=vPd1(:,3)==0;

if isempty(RC.Cr2)==1
    dPc=zeros(size(RC.Cr2,1),1);
    vPc1=zeros(size(RC.Cr2,1),4);
else
    dPc(:,1)=Pj(na+RC.Cr2)-Pj(na+RC.Cc2);
    DZ=dZ(2,:);
    dZW=DZ{1};  dZO=DZ{2};  dZG=DZ{3};
    vPc1(:,1)=dPc+dZW(RC.Cr2+(RC.Cc2-1)*nc)>=0;
    vPc1(:,2)=vPc1(:,1)==0;
    vPc1(:,3)=dPc+dZO(RC.Cr2+(RC.Cc2-1)*nc)>=0;
    vPc1(:,4)=vPc1(:,3)==0;
end

if isempty(RC.Gr2)==1
    dPg=zeros(size(RC.Gr2,1),1);
    vPg1=zeros(size(RC.Gr2,1),4);
else
    dPg(:,1)=Pj(na+nc+RC.Gr2)-Pj(na+nc+RC.Gc2);
    DZ=dZ(3,:);
    dZW=DZ{1};  dZO=DZ{2};  dZG=DZ{3};
    vPg1(:,1)=dPg+dZW(RC.Gr2+(RC.Gc2-1)*ng)>=0;
    vPg1(:,2)=vPg1(:,1)==0;
    vPg1(:,3)=dPg+dZO(RC.Gr2+(RC.Gc2-1)*ng)>=0;
    vPg1(:,4)=vPg1(:,3)==0;
end