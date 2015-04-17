function Q1=QBild(W1,W6,W7,Pi,uf,Won,dt,Pw,WonM,nw,wn)

dP=(Pi(Won)-Pw)*dt;
ql=W1.*dP;
qw=W6.*dP;
qp=W7.*dP;
Q=zeros(size(uf,1),5);

Q(uf==-1,1)=ql(uf==-1); 
Q(uf==1,2)=ql(uf==1); 
Q(uf==1,3)=ql(uf==1)-qw(uf==1); 
Q(uf==1,4)=0;%qg(uf==-1); 
Q(uf==1,5)=qp(uf==1); 
 
Q1=zeros(nw,5);
if numel(W1)~=0
    for i=1:5
        for j=1:size(wn,1)
        Q1(wn(j),i)=sum(Q(wn(j)==WonM,i));
        end
    end
end
