% √отовим координаты скважин и контуры
load('WellCoord')

WCoord2(:,1)=N2C(WCoord2(:,1));
WCoord2(:,1)=PrepText(WCoord2(:,1));
WCoord2(:,1)=PoiskTire(WCoord2(:,1));

de=zeros(size(wn));
WXY1=zeros(size(wn,1),2);
WXY2=zeros(size(wn,1),2);
WXY3=zeros(size(wn,1),2);


for i=1:size(wn,1)
    im=strcmp(wn(i),WCoord2(:,1));
    r=find(im);
    if isempty(r)==0
        de(i)=r(1);
        WXY2(i,1)=WCoord2{de(i),5};
        WXY2(i,2)=WCoord2{de(i),6};
    end;
end;


WCoord1(:,1)=N2C(WCoord1(:,1));
WCoord1(:,1)=PrepText(WCoord1(:,1));
WCoord1(:,1)=PoiskTire(WCoord1(:,1));
de1=zeros(size(wn));

for i=1:size(wn,1)
    im=strcmp(wn(i),WCoord1(:,1));
    r=find(im);
    if isempty(r)==0
        de1(i)=r(1);
        WXY1(i,1)=WCoord1{de1(i),3};
        WXY1(i,2)=WCoord1{de1(i),4};
    end;
end;

WCoord3(:,1)=N2C(WCoord3(:,1));
WCoord3(:,1)=PrepText(WCoord3(:,1));
WCoord3(:,1)=PoiskTire(WCoord3(:,1));
de2=zeros(size(wn));

for i=1:size(wn,1)
    im=strcmp(wn(i),WCoord3(:,1));
    r=find(im);
    if isempty(r)==0
        de2(i)=r(1);
        WXY3(i,1)=WCoord3{de2(i),2};
        WXY3(i,2)=WCoord3{de2(i),3};
    end;
end;