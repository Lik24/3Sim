function [ crack2 ] = Fracture( c,dl )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

for i=1:size(c,2)

    s=size(c{i},1);
    for j=1:s-1
        
        x1=c{i}(j,1);
        x2=c{i}(j+1,1);
        y1=c{i}(j,2);
        y2=c{i}(j+1,2);
        n=sqrt((x2-x1).^2+(y2-y1).^2)/dl;
        n=fix(n)+1;
        %   n-1 ���������� �������������� �����
        if n>=2   % ���� n>=2 �� ��������� ����, ����� ���������� ���� ��������� �����
            if (x1~=x2)*(y1~=y2)==1 % �� ������������ � �� �������������� �����
                k=(y1-y2)/(x1-x2);
                b=y1-x1*(y1-y2)/(x1-x2);
                for jj=1:n-1
                    a{i,j}(1,1)=x1;
                    a{i,j}(1,2)=y1;
                    a{i,j}(jj+1,1)=x1+jj*(x2-x1)/n;
                    a{i,j}(jj+1,2)=k*(x1+jj*(x2-x1)/n)+b;
                end;
            end;   % �� ������ � �� ����� �����
            
            
            if x1==x2 % ������������ �����
                for jj=1:n-1
                    a{i,j}(1,1)=x1;
                    a{i,j}(1,2)=y1;
                    a{i,j}(jj+1,1)=x1;
                    a{i,j}(jj+1,2)=y1+jj*(y2-y1)/n;
                end;
            end;   % ������������ �����
            
            if y1==y2 % �������������� �����
                for jj=1:n-1
                    a{i,j}(1,1)=x1;
                    a{i,j}(1,2)=y1;
                    a{i,j}(jj+1,1)=x1+jj*(x2-x1)/n;
                    a{i,j}(jj+1,2)=y1;
                end;
            end;   % �������������� �����
        else  a{i,j}(1,1)=x1; a{i,j}(1,2)=y1;
        end; %����� �����, ���� n>=2
        
    end;
    
end;

for i=1:size(a,1)
    zzz=[];
    for j=1:size(a,2)  zzz=[zzz;a{i,j}]; end;
    sz=size(zzz,1);
    sc=size(c{i},1);
    zzz(sz+1,1)=c{i}(sc,1);
    zzz(sz+1,2)=c{i}(sc,2);
    crack2{i,1}=zzz;
end;



