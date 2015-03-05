function [ pm2 ] = Purgatory( p,crack2,dl,WXY,rad,XY_GY)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% Первая часть чистит вокруг трещины
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
pm=[z,y];

% plot(pm(:,1),pm(:,2),'.')
[IN,ON]=inpolygon(pm(:,1),pm(:,2),XY_GY(:,1),XY_GY(:,2));

%Вторая часть чистит вокруг скважин

ww=zeros(size(pm,1),1);

%    for j =1:size(WXY,1)
       % xx=WXY(j,1:2);
        [k,d]=dsearchn(WXY,pm);
        ww(d<=rad)=1;
        ww(d==0)=0;
        ww(ON==1)=0;
        
zz(:,1)=pm(:,1);
yy(:,1)=pm(:,2);
zz(ww==1)=[];
yy(ww==1)=[];
pm2=[zz,yy];


end

