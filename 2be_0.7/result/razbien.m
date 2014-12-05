X=[0,2,0,2,9,10,9,10]';
Y=[0,0,2,2,0,0,2,2]';
XY=[X,Y];
XY=rand(20,2)*10;
n=size(XY,1);

[AA]=MR_Prop(XY,1);
plot(XY(:,1),XY(:,2),'*')
[L,~,~,H]=Geome3(AA,XY,ones(n,1),ones(n,1));

A1=AA;
A1=A1.*(L<500);

A1=A1.*(A1>0);
A1=A1-sparse(1:n,1:n,sum(A1));
b0=zeros(n,1);
A1=A1-sparse(1:n,1:n,1);
x0=b0'/A1;

fl=1;
k=0;
gr=zeros(n,1);
while fl==1
    k=k+1;
    b=b0;
    st=find(gr==0);
    b(st(1))=1;
    x=b'/A1;
    r=find(x~=x0);
    gr(r)=k;
    fl=sum(gr==0)~=0;
end

gk=unique(gr);
for i=1:size(gk,1)
  xy=XY(gr==gk(i),:);
  if size(xy,1)>2
  a=convhull(xy(:,1),xy(:,2));
  else
  a=1:size(xy,1);
  end
  plot(xy(:,1),xy(:,2),'*',xy(a,1),xy(a,2),'--')
  hold on
end
