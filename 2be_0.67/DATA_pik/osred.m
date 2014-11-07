function Y=osred(x,y,X)

n2=size(X,1);
n1=91;
n1=floor(n1)+1;
y1=zeros(n2*n1,size(y,2));
y1(1:size(y,1),:)=y;

for i=1:size(y,2)
    Y1=reshape(y1(:,i),n1,n2);
   for j=1:size(Y1,2) 
    A=Y1(:,j);
    Y(j,i)=mean(A(A~=0));
   end
end;