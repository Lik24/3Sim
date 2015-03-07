function VZL_VORONOI(DATA,Sw,p,WXY,uf,Won3)
xy=DATA.XY;
BND=DATA.BND;
ka=DATA.ka;
uf=uf(Won3(:,1));

[V,C]=VoronoiLimit(xy(:,1),xy(:,2),BND); 
%  DT=delaunayTriangulation(xy(:,1),xy(:,2),BND);
%      IO = isInterior(DT);
%      FG=DT.ConnectivityList(IO,:);
%      TR=triangulation(FG,xy(:,1),xy(:,2));
%      
% [V,C] = voronoiDiagram(TR);
%[V,C]=voronoin(xy); 
Sw0=zeros(size(ka,1),size(Sw,2));
Sw1(p)=Sw;
Sw0(ka==1)=Sw1;
Sw0(ka==0)=nan;


figure1=figure(45);
dx=max(xy(:,1))-min(xy(:,1));
dy=max(xy(:,2))-min(xy(:,2));
asx=dx/dy;

axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[asx 1 1]);
set(axes1,'PlotBoxAspectRatio',[asx 1 1])
hold(axes1,'all');

% axes1 = axes('Parent',figure1,'YGrid','on','XGrid','on',...
%     'Position',[0.13 0.11 0.651609195402299 0.815],...
%     'FontSize',14);

for k=1:length(C)
    
    VertCell = V(C{k},:);
    if sum(VertCell(:,1)>21200)==0
    ZZ=0*ones(size(VertCell(:,1)));
    patch(VertCell(:,1),VertCell(:,2),ZZ,1-Sw0(k),'FaceColor','flat')
    end
end
%set(gca,'CLim',[0 1])
colormap jet
hold on

gWXY=WXY(Won3(:,2)==1,:);
gWon3=Won3(Won3(:,2)==1,:);
uW=unique(gWon3(:,1));
for i=1:size(uW,1)
    plot(gWXY(uW(i)==gWon3(:,1),1),gWXY(uW(i)==gWon3(:,1),2),'black','LineWidth',1.5)
end

[uW1,ai,b1]=unique(Won3(:,1));
WXY1=WXY(ai,:);
uf1=uf(ai);

plot(WXY1(uf1==-1,1),WXY1(uf1==-1,2),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',3,...
    'Marker','v',...
    'LineStyle','none',...
    'Color',[0 0 1]);
plot(WXY1(uf1==1,1),WXY1(uf1==1,2),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',3,...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[1 0 0]);

hold off