A=rand(120);
x=rand(120,1);

tic
for i=1:1000
b=A*x;
end
toc

[r,c]=find(A);
vrc=r+100*(c-1);
tic
for i=1:1000
v=A(vrc).*x(c);
b1=accumarray(r,v);
end
toc

sum(abs(b-b1))