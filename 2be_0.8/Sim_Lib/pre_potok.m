function [vPa,vPc1,vPc2,vPg1,vPg2]=pre_potok(P,RC)

na=RC.na;
nc=RC.nc;
ng=RC.ng;

r=RC.Arc(:,1);
c=RC.Arc(:,2);

Pa=P(1:na);
Pc=P(na+1:na+nc);
Pg=P(na+nc+1:na+nc+ng);

vPa(:,1)=Pa(r)>=Pa(c);
vPa(:,2)=vPa(:,1)==0;

vPc1=Pc(RC.Cr)>=Pc(RC.Cc);
vPc2=vPc1==0;

% size(Pg)
% Pg(RC.Gr)

vPg1=Pg(RC.Gr)>=Pg(RC.Gc);
vPg2=vPg1==0;
%vPg=[vPg1,vPg2];
