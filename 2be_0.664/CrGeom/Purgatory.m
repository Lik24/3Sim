function [ pm2 ] = Purgatory( p,crack2,dl,WXY,rad )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% ������ ����� ������ ������ �������
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
% lklk

%������ ����� ������ ������ �������

ww=zeros(size(pm,1),1);

    for j =1:size(WXY,1)
        xx=WXY(j,1:2);
        [k,d]=dsearchn(xx,pm);
        ww(d<=rad)=1;
        ww(d==0)=0;
    end;

zz(:,1)=pm(:,1);
yy(:,1)=pm(:,2);
zz(ww==1)=[];
yy(ww==1)=[];
pm2=[zz,yy];
end

