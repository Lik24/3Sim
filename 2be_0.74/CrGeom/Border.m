function [ num_GR, cor_GR ] = Border( GR,pm2 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
s=size(GR,1)-1; %количество сегментов
% num_GR{s,1}=[];
e=1;
for i=1:s
    on=zeros(size(pm2,1),1);
    x1=GR(i,1);
    x2=GR(i+1,1);
    y1=GR(i,2);
    y2=GR(i+1,2);
    if (x1~=x2).*(y1~=y2) % не вертикальная и не горизонтальная линия
          k=(y1-y2)/(x1-x2);
          b=y1-x1*(y1-y2)/(x1-x2);
          for j=1:size(pm2,1)
              if (pm2(j,1)>=min([x1,x2])).*(pm2(j,1)<=max([x1,x2])).*(pm2(j,2)<=k*pm2(j,1)+b+e).*(pm2(j,2)>=k*pm2(j,1)+b-e) on(j,1)=1; end;
          end;
    end;
    if x1==x2 %Вертикальная линия
           for j=1:size(pm2,1)
               if (pm2(j,1)==x1).*(pm2(j,2)>=min([y1,y2])).*(pm2(j,2)<=max([y1,y2])) on(j,1)=1; end;
           end;
    end;
    
        if y1==y2 %Горизонтальная линия
           for j=1:size(pm2,1)
               if (pm2(j,2)==y1).*(pm2(j,1)>=min([x1,x2])).*(pm2(j,1)<=max([x1,x2])) on(j,1)=1; end;
           end;
    end;
    num_GR{i,1}=find(on);
end;
x_p=pm2(:,1);
y_p=pm2(:,2);

for i=1:s
   a=x_p(num_GR{i,1});
   b=y_p(num_GR{i,1});
   c(:,1)=a;
   c(:,2)=b;
   cor_GR{i,1}=c;
   clear a b c;
end;
end

