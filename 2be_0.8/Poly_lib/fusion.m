function wmu=fusion(Cpe,mup)
x1=mup(1,1);x2=mup(2,1);y1=mup(1,2);y2=mup(2,2);
% mu=[1,300,0.1,1];  
k=(y1-y2)/(x1-x2);
b=(y2*x1-y1*x2)/(x1-x2);


wmu=k*Cpe+b;

% wmu=((1-Cpe)./mu(1)+Cpe./mu(4));

end

