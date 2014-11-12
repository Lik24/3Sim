function VZL_VORONOI(xy,Sw,p,WXY,uf,BND)

[V,C]=VoronoiLimit(xy(:,1),xy(:,2),BND); 
%  DT=delaunayTriangulation(xy(:,1),xy(:,2),BND);
%      IO = isInterior(DT);
%      FG=DT.ConnectivityList(IO,:);
%      TR=triangulation(FG,xy(:,1),xy(:,2));
%      
% [V,C] = voronoiDiagram(TR);
%[V,C]=voronoin(xy); 

%Sw(p)=Sw;

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
    ZZ=0*ones(size(VertCell(:,1)));
    patch(VertCell(:,1),VertCell(:,2),ZZ,1-Sw(k),'FaceColor','flat')
    
end
%set(gca,'CLim',[0 1])
colormap jet
hold on

plot(WXY(uf==-1,1),WXY(uf==-1,2),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',4,...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
plot(WXY(uf==1,1),WXY(uf==1,2),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',4,...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[1 0 0]);
hold off