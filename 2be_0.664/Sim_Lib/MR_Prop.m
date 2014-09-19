function [A]=MR_Prop(WXY,Nl)
ns=size(WXY,1);

A=sparse(ns,ns);
cl=zeros(ns,1);

if ns==1
    Ri=inf;
elseif ns==2
    Ri=((WXY(1,1)-WXY(2,1)).^2+(WXY(1,2)-WXY(2,2)).^2).^0.5;
else
     DT=delaunayTriangulation(WXY(:,1),WXY(:,2));
%     triplot(DT)
%     kljgjh
    v2=1:ns;
    for i=1:ns
%         for j=1:ns
%             cl(j)=isConnected(DT,i,j);
%         end;
        v1=i*ones(ns,1);
        cl=isConnected(DT,v1,v2');
        A(i,cl~=0)=1;
    end;
end;

%A=repmat(A,Nl,Nl);
Ac=cell(Nl*Nl,1);
Ac(:)={sparse(ns,ns)};
Ac(1:Nl+1:end)={A};
Ac=reshape(Ac,Nl,Nl);
A=cell2mat(Ac);

v1=ns*ns*Nl+1:ns*Nl+1:ns*ns*Nl*Nl;
v2=ns+1:ns*Nl+1:ns*ns*Nl*Nl-Nl*ns*ns;

A2=sparse([v1,v2],ones(size([v1,v2])),ones(size([v1,v2])),ns*Nl*ns*Nl,1);
A2=reshape(A2,ns*Nl,ns*Nl);
A=A+A2-speye(size(A));