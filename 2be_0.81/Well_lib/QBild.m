function [Q1]=QBild(QQ,QQwom,uf,Won,dt,WonM,nw,W1)

Q=zeros(size(uf,1),5);

Q(uf==-1,1)=-QQ(uf==-1)*dt; 
Q(uf==1,2)=-QQ(uf==1)*dt; 
Q(uf==1,3)=-QQwom(uf==1,2)*dt; 
Q(uf==1,4)=-(QQ(uf==1,1) - QQwom(uf==1,1) - QQwom(uf==1,2))*dt; 
Q(uf==1,5)=0;%dt*qp(uf==1); 
 
wn=unique(WonM);
Q1=zeros(nw,5);
if numel(W1)~=0
    for i=1:5
        for j=1:size(wn,1)
        Q1(wn(j),i)=sum(Q(wn(j)==WonM,i));
        end
    end
end
