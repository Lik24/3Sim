function [vPa1,vPc1,vPg1,vPd1,dPa,dPc,dPg,dPd]=pre_potok_2(P,Pj,RC,rc_in_h,rc_in_hd,Na,Nd,na,nd,dZ,DZA)

%na=RC.na;
nc=RC.nc;
ng=RC.ng;

dPa=P(rc_in_h(:,2))-P(rc_in_h(:,1));
vPa1(:,1)=dPa+DZA(1).A(:,1)>=0;
vPa1(:,2)=vPa1(:,1)==0;

vPa1(:,3)=dPa+DZA(1).A(:,2)>=0;
vPa1(:,4)=vPa1(:,3)==0;

dPd=P(Na+nc+ng+rc_in_hd(:,2))-P(Na+nc+ng+rc_in_hd(:,1));
vPd1(:,1)=dPd+DZA(4).A(:,1)>=0;
vPd1(:,2)=vPd1(:,1)==0;

vPd1(:,3)=dPd+DZA(4).A(:,2)>=0;
vPd1(:,4)=vPd1(:,3)==0;

if isempty(RC.Cr2)==1
    dPc=zeros(size(RC.Cr2,1),1);
    vPc1=zeros(size(RC.Cr2,1),4);
else
    dPc(:,1)=Pj(na+RC.Cr2)-Pj(na+RC.Cc2);
    vPc1(:,1)=dPc+DZA(2).A(:,1)>=0;
    vPc1(:,2)=vPc1(:,1)==0;
    vPc1(:,3)=dPc+DZA(2).A(:,2)>=0;
    vPc1(:,4)=vPc1(:,3)==0;
end

if isempty(RC.Gr2)==1
    dPg=zeros(size(RC.Gr2,1),1);
    vPg1=zeros(size(RC.Gr2,1),4);
else
    dPg(:,1)=Pj(na+nc+RC.Gr2)-Pj(na+nc+RC.Gc2);
    vPg1(:,1)=dPg+DZA(3).A(:,1)>=0;
    vPg1(:,2)=vPg1(:,1)==0;
    vPg1(:,3)=dPg+DZA(3).A(:,2)>=0;
    vPg1(:,4)=vPg1(:,3)==0;
end