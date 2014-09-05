function [vPa1,vPc1,vPc2,vPg1,vPg2]=pre_potok_2(P,Pj,RC,rc_in_h,Na)

%na=RC.na;
nc=RC.nc;
%ng=RC.ng;

% r=RC.Arc(:,1);
% c=RC.Arc(:,2);
% 
%Pa=P(1:na);
% Pc=P(na+1:na+nc);
% Pg=P(na+nc+1:na+nc+ng);

vPa1(:,1)=P(rc_in_h(:,2))>=P(rc_in_h(:,1));
vPa1(:,2)=vPa1(:,1)==0;

vPc1=Pj(Na+RC.Cr2)>=Pj(Na+RC.Cc2);
vPc2=vPc1==0;

% size(Pg)
% Pg(RC.Gr)

vPg1=Pj(Na+nc+RC.Gr2)>=Pj(Na+nc+RC.Gc2);
vPg2=vPg1==0;
%vPg=[vPg1,vPg2];
