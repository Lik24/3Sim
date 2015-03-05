function [GY,WXY,H,Hk,K,Mp,Sw,Z1,Z3,GY_subl,Pw,Qw,Uf,Cpw,PwQC_bnd,Won3]=read_mr_prop
ff=pwd;
cd('DATA_In');
D=dir;
cd(ff);
    
for i=1:size(D,1)
    TT(i)={D(i).name};
end;

    mark(1)={'Bond_coord'};  mark(2)={'GY_sub'};     mark(3)={'Ha'};
    mark(4)={'Hk'}; mark(5)={'K'};   mark(6)={'Mp'};
    mark(7)={'Sw'};  mark(8)={'Well_coord2'};   mark(9)={'Z1'}; 
    mark(10)={'Z3'}; 
    mark(11)={'Pw'}; mark(12)={'Qw'}; mark(13)={'well_type'}; mark(14)={'Cpw'};
    
    for i=1:14
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
            WXY=table2array(T(:,2:3));
            A2=table2array(T(:,4:5));
            WN=table2array(T(:,1));
            Prd=ones(size(WN));
            
            v1=isnan(sum(A2,2))==0;
            if sum(v1)>0
                A2=A2(v1,:);
                wn=WN(v1,:);
                A1=table2array(T(:,2:3));
                WXY=[A1;A2];
                WN=[WN;wn];
                Prd=[Prd;2*ones(size(wn))];
            end
            [B,I]=sort(WN);
            WXY=WXY(I,:);
            Prd=Prd(I);
            WN(:,1)=B;
            uB=unique(B);
            g_flag=zeros(size(uB));

            for j=1:size(uB,1)
                if sum(WN(:,1)==uB(j))>1
                    g_flag(j)=1;
                end;
            end

            for j=1:size(uB,1)
                WN(WN(:,1)==uB(j),2)=g_flag(j);
            end
            WN(:,3)=Prd;
            A(i)={WXY};
            
        else
            T=readtable(TT{N(i)},'FileType','text','Delimiter','tab');
            A(i)={table2array(T)};
        end;
    end
    
    GY=A{1};
    WXY=A{8};
    Won3=WN;
    GY_subl=A{2};
    H=ChekFull(WXY,A{3});
    Hk=ChekFull(WXY,A{4});
    K=ChekFull(WXY,A{5});
    Mp=ChekFull(WXY,A{6});
    Sw=ChekFull(WXY,A{7});
    Z1=ChekFull(WXY,A{9});
    Z3=ChekFull(WXY,A{10});
    Pw1=A{11};Pw=time_data(Pw1);
    Qw1=A{12};Qw=time_data(Qw1);
    Uf1=A{13}; Uf=time_data(Uf1);
    Cpw1=A{14}; Cpw=time_data(Cpw1);
    
    B=importdata('Chek_bond.txt'); PwQC_bnd=B.data;
    

end
function A=ChekFull(WXY,B)
nb=size(B);
n=size(WXY,1);
if nb(1)<n
  for j=1:nb(2)
    A(:,j)=B(1,j)*ones(n,1);
  end;
else
    A=B;
end
end

function A=time_data(B)
N=size(B,1);
A=[];
for i=1:N
    A1=repmat(B(i,2:end)',1,B(i,1));
    A=[A,A1];
end;

end