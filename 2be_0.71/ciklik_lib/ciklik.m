function [qe]=ciklik(T,w)
global TT
n=size(T,1);
pz=10;
C=1e-4;
bet=0.01;
h1=ones(n,1);
h2=2*ones(n,1);
k1=ones(n,1);
k2=10*ones(n,1);

psi=bet*abs(k1-k2).*h1.*h2./(k1.*h1+k2.*h2);
B=ones(n,1);
A=ones(n,1);

vB=sparse(1:n,1:n,bet*w);
vA=sparse(1:n,1:n,-bet*w);
TT=T;
%x0=ones(2*n,1);

x=fsolve(@mufun,[bet*w*A;-bet*w*B]);
for i=1:10
vb=bet*w*B;
A=vb/T;
va=-bet*w*A;
B=va/T;
as(i,1)=sum(A);
as(i,2)=sum(B);
end

qe=psi*pz*C*w*((A.^2+B.^2)/2).^0.5;
qe
end
function F=mufun(x)
global TT
n=size(x,1);
F=zeros(n,1);
F(1:n/2)=TT*x(1:n/2);
F(n/2+1:end)=TT*x(n/2+1:end);
end