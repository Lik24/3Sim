function Q=Q2Sut(Qm,Qc,Qg,dt1,Ta)

for i=1:5
    qcs=size(Qc(:,i,:));
    QM(:,:)=Qm(:,i,:);
    QC(1:qcs(1),1:qcs(3))=Qc(:,i,:);
    
    qcs=size(Qg(:,i,:));
    QG(1:qcs(1),1:qcs(3))=Qg(:,i,:);

    Q(:,i,:)=[QM;QC;QG];
end;
size(Q)
%   Qo(:,1)=sum(Q(:,3,:));
%   Ql(:,1)=sum(Q(:,2,:));
%   Qz(:,1)=sum(Q(:,1,:));

for i=1:size(dt1,2)
 sQz(i,:)=sum(Q(:,1,1:i),3);
 sQl(i,:)=sum(Q(:,2,1:i),3);
 sQo(i,:)=sum(Q(:,3,1:i),3);
 sQg(i,:)=sum(Q(:,4,1:i),3);
 sQp(i,:)=sum(Q(:,5,1:i),3);
 st(i,1)=sum(dt1(1:i));
end;

sQz_d=interp1(st,sQz,1:Ta);
sQl_d=interp1(st,sQl,1:Ta);
sQo_d=interp1(st,sQo,1:Ta);
sQg_d=interp1(st,sQg,1:Ta);
sQp_d=interp1(st,sQp,1:Ta);


 Qz(1,:)=sQz_d(1,:);
 Ql(1,:)=sQl_d(1,:);
 Qo(1,:)=sQo_d(1,:);
 Qgaz(1,:)=sQg_d(1,:);
 Qp(1,:)=sQp_d(1,:);
 
 Qz(2:Ta,:)=sQz_d(2:Ta,:)-sQz_d(1:Ta-1,:);
 Ql(2:Ta,:)=sQl_d(2:Ta,:)-sQl_d(1:Ta-1,:);
 Qo(2:Ta,:)=sQo_d(2:Ta,:)-sQo_d(1:Ta-1,:);
 Qgaz(2:Ta,:)=sQg_d(2:Ta,:)-sQg_d(1:Ta-1,:);
 Qp(2:Ta,:)=sQp_d(2:Ta,:)-sQp_d(1:Ta-1,:);

Q=[];
Q(:,1,:)=Qz';
Q(:,2,:)=Ql';
Q(:,3,:)=Qo';
Q(:,4,:)=Qgaz';
Q(:,5,:)=Qp';