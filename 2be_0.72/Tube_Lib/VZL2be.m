function [PXY]=VZL2be(XY,XY2,r2,c2,p,DT,I1,I2,nt1,fl1,r,c)
figure(1), subplot(2,3,1), triplot(DT);

pX=[XY(r,1),XY(c,1)];
pY=[XY(r,2),XY(c,2)];

figure(1);  subplot(2,3,2),
for i=1:size(r,1)
    plot(pX(i,:),pY(i,:))
    hold on
end;
PXY(1,1)={pX};
PXY(2,1)={pY};

pX=[XY2(r2,1),XY2(c2,1)];
pY=[XY2(r2,2),XY2(c2,2)];

figure(1),  subplot(2,3,2),
for i=1:size(r2,1)
    plot(pX(i,:),pY(i,:),'g')
    hold on
end;
plot(XY(nt1,1),XY(nt1,2),'r*')

xy=XY(fl1==1,:);
plot(xy(I1,1),xy(I1,2),'mo');
plot(xy(I2,1),xy(I2,2),'mo');


cp=p/max(p);

 subplot(2,3,3)
% for i=1:n
% plot(xy(i,1),xy(i,2),'o','MarkerFaceColor',[cp(i),cp(i),cp(i)])
% hold on
% end;
% tri = delaunay(xy(:,1),xy(:,2));
% trisurf(tri,xy(:,1),xy(:,2),p)
xg=XY(:,1);
xg(fl1==1)=[];
yg=XY(:,2);
yg(fl1==1)=[];
pv=zeros(1,size(yg,1));
pv(:)=NaN;
[X,Y]=meshgrid(0:1:1000,0:1:1000);
F=scatteredInterpolant([xy(:,1);xg],[xy(:,2);yg],[p,pv]','linear','none');
P=F(X,Y);
contourf(X,Y,P)
hold on
[C,I1]=min(XY(fl1==1,1)+XY(fl1==1,2));
[C,I2]=max(XY(fl1==1,1)+XY(fl1==1,2));
xy=XY(fl1==1,:);
plot(xy(I1,1),xy(I1,2),'mo');
plot(xy(I2,1),xy(I2,2),'mo');
for i=1:size(r2,1)
    plot(pX(i,:),pY(i,:),'g')
    hold on
end;
plot(XY(nt1,1),XY(nt1,2),'r*')

 subplot(2,3,4)
for i=1:size(r2,1)
    plot(pX(i,:),pY(i,:),'Color',[cp(r2(i)),cp(r2(i)),0.5],'LineWidth',4)
    hold on
end;
