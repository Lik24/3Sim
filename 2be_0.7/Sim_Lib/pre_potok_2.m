function [vPa1,vPc1,vPc2,vPg1,vPg2,vPd1,dPa,dPc,dPg,dPd]=pre_potok_2(P,Pj,RC,rc_in_h,rc_in_hd,Na,Nd,na,nd)

%na=RC.na;
nc=RC.nc;
ng=RC.ng;

% r=RC.Arc(:,1);
% c=RC.Arc(:,2);
% 
%Pa=P(1:na);
% Pc=P(na+1:na+nc);
% Pg=P(na+nc+1:na+nc+ng);

dPa=P(rc_in_h(:,2))-P(rc_in_h(:,1));

vPa1(:,1)=dPa>=0;
vPa1(:,2)=vPa1(:,1)==0;

dPd=P(Na+nc+ng+rc_in_hd(:,2))-P(Na+nc+ng+rc_in_hd(:,1));
vPd1(:,1)=dPd>=0;
vPd1(:,2)=vPd1(:,1)==0;

dPc=Pj(na+RC.Cr2)-Pj(na+RC.Cc2);
vPc1=dPc>=0;
vPc2=vPc1==0;

dPg=Pj(na+nc+RC.Gr2)-Pj(na+nc+RC.Gc2);
vPg1=dPg>=0;
vPg2=vPg1==0;
