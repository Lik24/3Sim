function [vPa1,vPc1,vPc2,vPg1,vPg2,dPa,dPc,dPg]=pre_potok_2(P,Pj,RC,rc_in_h,Na)

%na=RC.na;
nc=RC.nc;
%ng=RC.ng;

% r=RC.Arc(:,1);
% c=RC.Arc(:,2);
% 
%Pa=P(1:na);
% Pc=P(na+1:na+nc);
% Pg=P(na+nc+1:na+nc+ng);

dPa=P(rc_in_h(:,2))-P(rc_in_h(:,1));

vPa1(:,1)=dPa>=0;
vPa1(:,2)=vPa1(:,1)==0;

dPc=Pj(Na+RC.Cr2)-Pj(Na+RC.Cc2);
vPc1=dPc>=0;
vPc2=vPc1==0;

dPg=Pj(Na+nc+RC.Gr2)-Pj(Na+nc+RC.Gc2);
vPg1=dPg>=0;
vPg2=vPg1==0;
%vPg=[vPg1,vPg2];
