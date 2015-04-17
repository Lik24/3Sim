function [ pm2 ] = Purgatory( p,crack2,dl,WXY,rad,fl1,fl2)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% Первая часть чистит вокруг трещины
if fl1==1
w=zeros(size(p,1),1);
for i=1:size(crack2,1)
    for j =1:size(crack2{i,1},1)
        x=crack2{i,1}(j,1:2);
        [k,d]=dsearchn(x,p);
        w(d<=dl)=1;
    end;
end;
z(:,1)=p(:,1);
y(:,1)=p(:,2);
z(w==1)=[];
y(w==1)=[];
p=[z,y];
end
% plot(pm(:,1),pm(:,2),'.')
% lklk

%Вторая часть чистит вокруг скважин

if fl2==1
ww=zeros(size(p,1),1);

%    for j =1:size(WXY,1)
       % xx=WXY(j,1:2);
        [k,d]=dsearchn(WXY,p);
        ww(d<=rad)=1;
        ww(d==0)=0;
%    end;

zz(:,1)=p(:,1);
yy(:,1)=p(:,2);
zz(ww==1)=[];
yy(ww==1)=[];
pm2=[zz,yy];
end

end

