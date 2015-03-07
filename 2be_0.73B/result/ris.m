%load('matlab16.2.mat');

for i=1:size(SD1,2)
  CD=SD{i};
  Q=CD{1};
  Qo(:,i)=sum(Q(:,3,:));
  Ql(:,i)=sum(Q(:,2,:));
  Qz(:,i)=sum(Q(:,1,:));
end;

for i=1:size(Qo,1)
 sQo(i,:)=sum(Qo(1:i,:),1);
 sQl(i,:)=sum(Ql(1:i,:),1);
end;

c=1-Qo./Ql;