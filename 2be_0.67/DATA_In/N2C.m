function A=N2C(A)
for i=1:size(A,1)
    if ischar(A(i))==0
        A(i)={num2str(A{i})};
    end;
end;
end




