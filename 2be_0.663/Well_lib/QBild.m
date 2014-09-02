function Q=QBild(W1,W6,W7,Pi,uf,Won,dt,Pw)

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
 

