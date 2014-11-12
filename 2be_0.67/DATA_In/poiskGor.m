% 
% for l=1:3
% for i=1:1098
%     if isempty(H{i,l})==0
%         h(i,l)=H{i,l};
%     else
%         h(i,l)=0;
%     end
% end
% end


f0=[]; % поиск горизонтальных скважин
for i=1:size(WN,1)
    if ~isempty(findstr('H',WN{i,1}))
        f0=[f0;i];
    end;
    
    if ~isempty(findstr('S',WN{i,1}))
        f0=[f0;i];
    end;
end;


[A]=MR_Prop(WXY,1);
H1=H;
H1k=Hk;
f1=[];

for i=1:size(f0,1)
    A1=A(:,f0(i));
    A1(f0)=0;
    c=find(A1);
    if isempty(c)==0
        for l=1:3
            H1(f0(i),l)=mean(H(c,l));
            H1k(f0(i),l)=mean(Hk(c,l));
        end;
    else
        f1=[f1,f0(i)];
    end
end;
% 
for i=1:size(f1,1)
    A1=A(:,f1(i));
    A1(f1)=0;
    c=find(A1);

    for l=1:3
      H1(f1(i),l)=mean(H1(c,l));
      H1k(f1(i),l)=mean(H1k(c,l));
    end;
end;