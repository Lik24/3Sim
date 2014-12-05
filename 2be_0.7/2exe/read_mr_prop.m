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
    H=ChekFull(WXY,A{3});
    Hk=ChekFull(WXY,A{4});
    K=ChekFull(WXY,A{5});
    Mp=ChekFull(WXY,A{6});
    Sw=ChekFull(WXY,A{7});
    Z1=ChekFull(WXY,A{9});
    Z3=ChekFull(WXY,A{10});
end
function A=ChekFull(WXY,B)
nb=size(B);
n=size(WXY,1);
if nb(1)<n
  for j=1:nb(2)
    A(:,j)=B(1,j)*ones(n,1);
  end;
else
    A=B{:};
end
end