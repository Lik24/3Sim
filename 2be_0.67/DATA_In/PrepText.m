function A=PrepText(B)
    A=cell(size(B(:)));
    for i=1:size(B(:),1)
     A(i)={deblank(B{i})};
     A(i)={upper(A{i})};
    end;
    A=reshape(A,size(B));
end