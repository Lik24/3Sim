function [GY,WXY,H,Hk,K,Mp,Sw,Z1,Z3,GY_subl,Pw,Qw,Uf,Cpw,PwQC_bnd,Won3]=read_mr_prop_MF1

    mark(1)={'Bond_coord'};  mark(2)={'GY_sub'};     mark(3)={'Ha'};
    mark(4)={'Hk'}; mark(5)={'K'};   mark(6)={'Mp'};
    mark(7)={'Sw'};  mark(8)={'Well_coord2'};   mark(9)={'Z1'}; 
    mark(10)={'Z3'}; 
    mark(11)={'Pw'}; mark(12)={'Qw'}; mark(13)={'well_type'}; mark(14)={'Cpw'};
%%  Координаты границы участка  
A(1)={[0,    0;
      0,    500;
      500, 500;
      500, 0;
      0,    0]};
 %%  Координаты границы слоёв 
 Nl=2;
 a2=A{1};
 a2=[1,-1;a2];
 a2=repmat(a2,Nl,1);
 a2(1:6:end,1)=1:Nl;
 A(2)={a2};
%  A(2)={[1,    -1;
%      0,    0;
%      0,    500;
%      500, 500;
%      500, 0;
%      2,    -1;
%      0,	0;
%      0,	500;
%      500,	500;
%      500,	0;
%      3,    -1;
%      0,	0;
%      0,	500;
%      500,	500;
%      500,	0;]};
%%  Значения параметров пласта по скважинам/слоям  
A(3)={10*ones(1,Nl)};        %  Толшина слоя  
A(4)={10*ones(1,Nl)};        %  Толшина коллектора  
A(5)={1000*ones(1,Nl)};  %  Проницаемость  
A(6)={0.30*ones(1,Nl)};  %  Пористость
A(7)={0.25*ones(1,Nl)};    %  Водонасыщенность
A(9)={repmat(800:2.5:800+2.5*(Nl-1),5,1)};      %  Глубина залегания
A(10)={repmat(800:2.5:800+2.5*(Nl-1),5,1)};      %  Глубина залегания
   
%% Значение для скважин
A(8)={[1, 0,   0,   nan,  nan;         %Координаты
       2, 0,   500, nan,  nan;
       3, 500, 500, nan,  nan;
       4, 500,   0, nan,  nan;
       5, 250,  250, nan, nan;
       6, 250,     0, nan, nan;
       7,   0,    250, nan, nan;
       8, 500,   250, nan, nan;
       9, 250,   500, nan, nan;]};

A(11)={[7300,20,20,400,20,20,20,20,20,20]};    %  Забойное давление
A(12)={[7300,0,0,0,0,-100*0,0,0,0,0]};           %  Дебит жидкости
A(13)={[7300,1,1,-1,1,1,1,1,1,1]};       %  Тип скважины
A(14)={[7300,0,0,0,0,0,0,0,0,0]};           %  Полимер
 
%% Запись
WXY=A{8}(:,2:3);
A2=A{8}(:,4:5);
WN=A{8}(:,1);
Prd=ones(size(WN));
v1=isnan(sum(A2,2))==0;
if sum(v1)>0
    A2=A2(v1,:);
    wn=WN(v1,:);
    A1=A{8}(:,2:3);
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
A(8)={WXY};
 
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
    
%% Технологические ограничения Pw_min	Pw_max	Qp_min	Qp_max	Qi_min	Qi_max	WCT	Qoil_min   
PwQC_bnd=[1	3000	0	200	-200	0	1.98	0;
          1	 70     1	200	-200	0	1.98	0;
          1	 70     1	200	-200	0	1.98	0;
          1	 70     1	200	-200	0	1.98	0;
          1	 70     1	200	-200	0	1.98	0;
          1	 70     1	200	-200	0	1.98	0;
          1	 70     1	200	-200	0	1.98	0;
          1	 70     1	200	-200	0	1.98	0;
          1	 70     1	200	-200	0	1.98	0];
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