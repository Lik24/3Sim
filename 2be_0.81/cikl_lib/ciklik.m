function [qe,A,B]=ciklik(flag,T,PR,Won,Uf,W1,na,b1gm,dv,Psi,MSw)

if (PR.zc(1)+PR.zc(2))==0
    PR.zc(1)=1e-6;
    PR.zc(2)=1e-5;
end
mC=MSw*PR.zc(1)+(1-MSw)*PR.zc(2); % средн€€ сжимаемость жидкости
 
if flag==1
    w=2*pi/PR.Tcikl;
    n=size(T,1);
   
    ba=zeros(n,1);
    bb=zeros(n,1);
    bb(Won(Uf(Won(:,3))==-1,1))=W1(Uf(Won(:,3))==-1,1)*1;
    bab=[ba;bb];
    
    vB=sparse(1:n,1:n,-mC.*w.*dv);
    vA=sparse(1:n,1:n,mC.*w.*dv);
    %TT=[T,vB;vA,T];
    TT=[T,vB;vA,T];
    
    x=bab'/TT;
    x=x';
    
    A=x(1:end/2);
    B=x(end/2+1:end);
    %qe=Pa*w*psi.*((A.^2+B.^2)/2).^0.5;
    qe=PR.P_Amp*w*Psi*((A.^2+B.^2)/2).^0.5;
    qe=qe.*mC;
else
    qe=zeros(na,1);
    A=zeros(na,1);
    B=zeros(na,1);
end