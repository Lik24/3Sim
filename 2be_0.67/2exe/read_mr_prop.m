function [GY,WXY,H,Hk,K,Mp,Sw,Z1,Z3,GY_subl]=read_mr_prop
ff=pwd;
cd('2exe');
D=dir;
cd(ff);
    
for i=1:size(D,1)
    TT(i)={D(i).name};
end;

    mark(1)={'Bond_coord'};  mark(2)={'GY_sub'};     mark(3)={'Ha'};
    mark(4)={'Hk'}; mark(5)={'K'};   mark(6)={'Mp'};
    mark(7)={'Sw'};  mark(8)={'Well_coord'};   mark(9)={'Z1'}; 
    mark(10)={'Z3'};
    
    for i=1:10
        for j=1:size(TT,2)
            k=findstr(mark{i},TT{j});
            if isempty(k)==0
                N(i)=j;
            end;
        end;
    end;

    for i=1:size(N,2)
        if i==8
            T=readtable(TT{N(i)},'Delimiter','tab');
            A(i)={table2array(T(:,2:3))};
        else
            T=readtable(TT{N(i)},'FileType','text','Delimiter','tab');
            A(i)={table2array(T)};
        end;
    end
    
    GY=A{1};
    WXY=A{8};
    GY_subl=A{2};
    H=A{3};
    Hk=A{4};
    K=A{5};
    Mp=A{6};
    Sw=A{7};
    Z1=A{9};
    Z3=A{10};
    