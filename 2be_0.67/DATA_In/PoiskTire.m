function A=PoiskTire(A)
for i=1:size(A,1)
A(i)=strrep(A(i),'_','-');
end;