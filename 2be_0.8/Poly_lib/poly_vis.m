function [TP,wmu]=poly_vis(Cpe,mup,fp,Kf,T,n,r,c,mu)
x1=mup(1,1);x2=mup(2,1);y1=mup(1,2);y2=mup(2,2);
k=(y1-y2)/(x1-x2);
b=(y2*x1-y1*x2)/(x1-x2);


if fp==1
    wmu=k*Cpe+b;
    Tp=T.*Kf.*Cpe./wmu;
    Tp1=sparse(r,c,Tp,n,n);
    TP=Tp1-sparse(1:n,1:n,sum(Tp1,2),n,n);
end;


if fp==0
    wmu=mu;
    TP=sparse(n,n);
end;


% Tp=T.*Kf.*Cpe./wmu;
%wmu=((1-Cpe)./mu(1)+Cpe./mu(4));

end

