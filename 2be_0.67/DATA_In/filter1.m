function [K,Sw,Mp,Hk,H]=filter1(K,Sw,Mp,Hk,H,WXY)

Ki=K~=0;
Swi=Sw~=0;
Mpi=Mp~=0;
Hki=Hk~=0;
Hi=H~=0;

V(:,:,1)=K;
V(:,:,2)=Sw;
V(:,:,3)=Mp;
V(:,:,4)=Hk;
V(:,:,5)=H;

%% Толшины
H(Hi==0)=Hk(Hi==0);

Si=(Ki+Swi+Mpi+Hki+Hi);

[r,c]=find((Si<5).*(Si~=0));

for i=1:size(r,1)
  v1(i,1)=sum(Ki(r(i),c(i)));
  v1(i,2)=sum(Swi(r(i),c(i)));
  v1(i,3)=sum(Mpi(r(i),c(i)));
  v1(i,4)=sum(Hki(r(i),c(i)));
  v1(i,5)=sum(Hi(r(i),c(i)));
end

[row,col]=find(v1==0);
for i=1:size(col,1)
     V(:,:,col(i))=Vostanovlenie(V(:,:,col(i)),WXY,r(row(i)),c(row(i)));
end;

K=V(:,:,1);
Sw=V(:,:,2);
Mp=V(:,:,3);
Hk=V(:,:,4);
H=V(:,:,5);

end

function V=Vostanovlenie(V,WXY,r,c)
cu=unique(c);
for l=cu'
    x=WXY(:,1);
    y=WXY(:,2);
    k=V(:,l);
    x(r(c==l))=[];
    y(r(c==l))=[];
    k(r(c==l))=[];
    x(k==0)=[];
    y(k==0)=[];
    k(k==0)=[];
    if size(x,1)>2
    F=scatteredInterpolant(x,y,k,'linear','nearest');
    V(r(c==l),l)=F(WXY(r(c==l),1),WXY(r(c==l),2));
    else
    V(r(c==l),l)=mean(k);    
    end
end;
end