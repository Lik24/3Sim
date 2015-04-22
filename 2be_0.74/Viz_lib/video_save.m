function video_save(DATA,XY,WXY,Z,Pi,Sw,Nl,pp,uf,SwC,NT,CR_GRUP,z,A2C,pc)
tic
writerObj = VideoWriter('SW_plot.avi','MPEG-4');
open(writerObj);
size(Pi)
size(XY)
xy=DATA.XY;
BND=DATA.BND;
ka=DATA.ka;
A2C=A2C(:,pc);

Sw0=zeros(size(ka,1),size(Sw,2));
Sw1=Sw(1:size(pp,2),:);
Sw0(ka==1,:)=Sw1;
Sw0(ka==0,:)=nan;
Sw0(pp,:)=Sw0(:,:);


P2=Pi(1:size(XY,1),:);
P3=Sw(1:size(XY,1),:);

XY(pp,1)=XY(:,1);
XY(pp,2)=XY(:,2);
t=size(P2,2);

% x=reshape(XY(:,1),size(Z,1)/Nl,Nl);
% y=reshape(XY(:,2),size(Z,1)/Nl,Nl);
x=XY(:,1);
y=XY(:,2);

mx(1)=min(XY(:,1));
mx(2)=max(XY(:,1));
my(1)=min(XY(:,2));
my(2)=max(XY(:,2));
[X,Y]=meshgrid(mx(1):10:mx(2),my(1):10:my(2));

flt=ones(2,t*10);
dy=365;
for t1=1:t*10
    if mod(t1,dy)<=0.5*dy
    flt(1,t1)=0;
    flt(2,t1)=1;
    else
    flt(1,t1)=1;
    flt(2,t1)=0;
    end;
end;
flt=flt(:,1:10:end);


for j=1:t
 [frame(j)]=gogi(P2(:,j),P3(:,j),pp,Z,Nl,xy(:,1),xy(:,2),WXY,X,Y,uf(:,j),SwC(:,j),NT,CR_GRUP,XY,BND,z,A2C,xy);
%frame = imdilate(frame, se);

end

for j=1:t
writeVideo(writerObj,frame(j));
end;
close(writerObj);
%set(gcf, 'Res', 'auto')
%print('-djpeg','-r300')
%saveas(h,strcat('pic/',num2str(ni),'/Карта_вар_',L,'.jpg')); 
toc
end

function [frame]=gogi(P2,P3,pp,WZ,Nl,x,y,WXY,X,Y,uf,SwC,NT,CR_GRUP,XY,BND,z,A2C,xy)
Pt(pp)=P2;
%p1=reshape(Pt,size(WZ,1)/Nl,Nl);
p1=Pt';
Pt(pp)=P3;
%p2=reshape(Pt,size(WZ,1)/Nl,Nl);
p2=Pt';

for i=1:Nl
F=scatteredInterpolant(x(:,i),y(:,i),p2(:,i),'linear','none');
B(:,:,i)=F(X,Y);
F=scatteredInterpolant(x(:,i),y(:,i),p1(:,i),'linear','none');
Bp(:,:,i)=F(X,Y);
end;

%h=figure(1);
% contourf(X,Y,mean(B,3),'LineStyle','none','LineColor',[0 0 0],...
%     'Fill','on');
h=plot_sw_p(X,Y,B,Bp,WXY(uf==-1,1),WXY(uf==-1,2),WXY(uf==1,1),WXY(uf==1,2),SwC,NT,Nl,CR_GRUP,XY,z,A2C,pp);
alpha 0.9
[V,C,I]=VoronoiLimit(xy(:,1),xy(:,2),BND); 
Sw0=P3(I,:);
for k=1:length(C)
    
    VertCell = V(C{k},:);
    ZZ=0*ones(size(VertCell(:,1)));
    subplot(1,1,1),patch(VertCell(:,1),VertCell(:,2),ZZ,Sw0(k,1),'FaceColor','none','EdgeColor',[0.6 0.6 0.6]);
    %subplot(1,2,2),patch(VertCell(:,1),VertCell(:,2),ZZ,Sw0(k,1),'FaceColor','none','EdgeColor',[0.6 0.6 0.6]);
end

%h=plot_sw(X,Y,B,WXY([1,4,7],1),WXY([1,4,7],2),WXY([3,6,9],1),WXY([3,6,9],2),PXY);
%h1=plot_p(X,Y,Bp,WXY([1,4,7],1),WXY([1,4,7],2),WXY([3,6,9],1),WXY([3,6,9],2),PXY);
% set(gca,'nextplot','replacechildren');
% set(gcf,'Renderer','zbuffer');
%set(gcf,'-r300');
frame = getframe(h);
close(h)
end