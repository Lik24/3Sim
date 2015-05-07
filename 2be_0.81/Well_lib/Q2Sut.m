function [Q,pw,ppl,dtz]=Q2Sut(Qm,Qc,Qg,Qd,Pw,PpW,dt1,Ta)

if Ta<10*365
    dtz=1;
elseif Ta<200*365
    dtz=365/12;
else
    dtz=365;
end;

for i=1:5
    qcs=size(Qc(:,i,:));
    QM(:,:)=Qm(:,i,:);
    QC(1:qcs(1),1:qcs(3))=Qc(:,i,:);
    
    qcs=size(Qg(:,i,:));
    QG(1:qcs(1),1:qcs(3))=Qg(:,i,:);

    qcs=size(Qd(:,i,:));
    QD(1:qcs(1),1:qcs(3))=Qd(:,i,:);
    
    Q(:,i,:)=[QM;QC;QG;QD];
end;

nt=size(dt1,2);

st(:,1)=cumsum(dt1);
sQ=cumsum(Q(:,:,1:nt),3);
sQz(:,:)=sQ(:,1,:);
sQl(:,:)=sQ(:,2,:);
sQo(:,:)=sQ(:,3,:);
sQg(:,:)=sQ(:,4,:);
sQp(:,:)=sQ(:,5,:);

sQz_d=interp1(st,sQz',dtz:dtz:Ta,'linear','extrap');
sQl_d=interp1(st,sQl',dtz:dtz:Ta,'linear','extrap');
sQo_d=interp1(st,sQo',dtz:dtz:Ta,'linear','extrap');
sQg_d=interp1(st,sQg',dtz:dtz:Ta,'linear','extrap');
sQp_d=interp1(st,sQp',dtz:dtz:Ta,'linear','extrap');

%sQz_d

 Qz(1,:)=sQz_d(1,:);
 Ql(1,:)=sQl_d(1,:);
 Qo(1,:)=sQo_d(1,:);
 Qgaz(1,:)=sQg_d(1,:);
 Qp(1,:)=sQp_d(1,:);
 
 Tend=size(sQz_d,1);
 
 Qz(2:Tend,:)=sQz_d(2:end,:)-sQz_d(1:end-1,:);
 Ql(2:Tend,:)=sQl_d(2:end,:)-sQl_d(1:end-1,:);
 Qo(2:Tend,:)=sQo_d(2:end,:)-sQo_d(1:end-1,:);
 Qgaz(2:Tend,:)=sQg_d(2:end,:)-sQg_d(1:end-1,:);
 Qp(2:Tend,:)=sQp_d(2:end,:)-sQp_d(1:end-1,:);


Q=[];
Q(:,1,:)=Qz';
Q(:,2,:)=Ql';
Q(:,3,:)=Qo';
Q(:,4,:)=Qgaz';
Q(:,5,:)=Qp';

ppl=interp1(st,PpW(:,1:end)',dtz:dtz:Ta,'linear','extrap')';
pw=interp1(st,Pw(:,1:end)',dtz:dtz:Ta,'linear','extrap')';