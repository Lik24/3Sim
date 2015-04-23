function [Q1]=QBild2(Wf,Won,Kfw,Kfo,mu,Uf,Pi,Pw,dt,WonM,nw,wn)

TW=Kfw(Won(:,1))/mu(1);
TO=Kfo(Won(:,1))/mu(2);

Tiw=(1-0)./mu(1);
Tip=0./mu(1);


W5=Wf.*TO.*(Uf==1);%
W6=Wf.*(TW.*(Uf==1)+(Tiw+Tip).*(Uf==-1));%./Bwo(Won,2)

dP=(Pi(Won(:,1))-Pw)*dt;
qo=W5.*dP;
qw=W6.*dP;

Q=zeros(size(Uf,1),5);

Q(Uf==-1,1)=qw(Uf==-1); 
Q(Uf==1,2)=qo(Uf==1)+qw(Uf==1); 
Q(Uf==1,3)=qo(Uf==1); 
Q(Uf==1,4)=0;%qg(uf==-1); 
%Q(uf==1,5)=qp(uf==1); 
 
Q1=zeros(nw,5);
if numel(W5)~=0
    for i=1:5
        for j=1:size(wn,1)
        Q1(wn(j),i)=sum(Q(wn(j)==WonM,i));
        end
    end
end
