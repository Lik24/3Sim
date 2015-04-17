function VZL2(DATA,WXY,P,Sw,T,Cp,Nl,pp,Q,SwC,CR_GRUP,pc,NT,XYgy,a0,pb,GYData,XYgy2,dtz,Won3,A2C,GY_subl)

XY=DATA.XY;
WZ=DATA.gZ(:);
ka=DATA.ka;
K=DATA.gKX;

a0=a0(Won3(:,1));

XY=repmat(XY,Nl,1);

P1=P(1:size(pp,2),:);
PC=P(size(pp,2)+1:size(pp,2)+size(pc,2),:);
T1=T(1:size(pp,2),:);
%K(pp)=K;
P1(pp,:)=P1;
Sw1(pp,:)=Sw;
Cp1(pp,:)=Cp;
T1(pp,:)=T1;
PC(pc,:)=PC;
SwC(pc,:)=SwC;
A2C(:,pc)=A2C;

%P1GY=P(size(pp,2)+1:end,:);
%P1GY(pb,:)=P1GY;
XY_GY2=GYData.XY;

XYgy=DATA.XY(DATA.BND(:,1),:);

P=zeros(size(ka,1),size(P,2));
T=zeros(size(ka,1),size(P,2));
Sw=zeros(size(ka,1),size(Sw,2));
Cp=zeros(size(ka,1),size(P,2));

P(ka==1,:)=P1;
T(ka==1,:)=T1;
Sw(ka==1,:)=Sw1;
Cp(ka==1,:)=Cp1;

P(ka==0,:)=nan;
T(ka==0,:)=nan;
Sw(ka==0,:)=nan;
Cp(ka==0,:)=nan;
K(ka==0)=nan;

nxy=size(WZ,1)/Nl;
%[A]=MR_Prop(XY,1);
%A=A.*(A>0);
for l=1:Nl
    
Z(nxy*(l-1)+1:nxy*l)=(-l+1)*20;
end;

x=reshape(XY(:,1),size(WZ,1)/Nl,Nl);
y=reshape(XY(:,2),size(WZ,1)/Nl,Nl);
% x1=reshape(XY_GY2(:,1),size(XY_GY2,1)/Nl,Nl);
% y1=reshape(XY_GY2(:,2),size(XY_GY2,1)/Nl,Nl);

z=reshape(Z,size(WZ,1)/Nl,Nl);
log10k=reshape((K/8.64),size(WZ,1)/Nl,Nl);
p=reshape(P(:,end),size(WZ,1)/Nl,Nl);
sw=reshape(Sw(:,end),size(WZ,1)/Nl,Nl);
cp=reshape(Cp(:,end),size(WZ,1)/Nl,Nl);
tt=reshape(T(:,end),size(WZ,1)/Nl,Nl);

%pgy=reshape(P1GY(:,end),size(XY_GY2,1)/Nl,Nl);

mx(1)=min(XY(:,1));
mx(2)=max(XY(:,1));
my(1)=min(XY(:,2));
my(2)=max(XY(:,2));
[X,Y]=meshgrid(mx(1):5:mx(2),my(1):5:my(2));

% mx(1)=min(XY_GY2(:,1));
% mx(2)=max(XY_GY2(:,1));
% my(1)=min(XY_GY2(:,2));
% my(2)=max(XY_GY2(:,2));
% [X1,Y1]=meshgrid(mx(1):5:mx(2),my(1):5:my(2));
%  
figure(98),s1=subplot(2,4,1);
% plot_fild(x1,y1,z,pgy,Nl,X1,Y1,WXY,'Пластовое давление',XYgy2,a0,'nearest') % 
% hold on
plot_fild(x,y,z,p,Nl,X,Y,WXY,GY_subl,'Пластовое давление',XYgy,a0,'nearest') % 
hold on
%ax=[XY_GY2(GYData.BND(:,1),1),XY_GY2(GYData.BND(:,2),1)];
%ay=[XY_GY2(GYData.BND(:,1),2),XY_GY2(GYData.BND(:,2),2)];
%plot(ax,ay,'k','LineWidth',2)
PCC=(abs(min(PC(:,end)))+PC(:,end))/max(abs(min(PC(:,end)))+PC(:,end));
if numel(PCC)~=0
plot_crack_color(Nl,NT,PCC,CR_GRUP,XY,z,A2C);
end
%set(s1,'CLim',[min([pgy;p]) max([pgy;p])])
hold off

figure(98),subplot(2,4,2);
plot_fild(x,y,z,sw,Nl,X,Y,WXY,GY_subl,'Водонасыщенность',XYgy,a0,'linear') % 
hold on
if numel(SwC)~=0
plot_crack_color(Nl,NT,SwC,CR_GRUP,XY,z,A2C);
end
hold off


figure(98),subplot(2,4,4);
plot_fild(x,y,z,cp,Nl,X,Y,WXY,GY_subl,'Концентрация',XYgy,a0,'nearest') % 

figure(98),subplot(2,4,3);
plot_fild(x,y,z,tt,Nl,X,Y,WXY,GY_subl,'Температура',XYgy,a0,'nearest') % 

figure(98),subplot(2,4,5);
plot_fild(x,y,z,log10k,Nl,X,Y,WXY,GY_subl,'Проницаемость',XYgy,a0,'nearest') % 


 subplot(2,4,8);
 T=(1:size(Q,3))*dtz;
% 
Qz(:,:)=sum(Q(:,1,:))/dtz;
Qd(:,:)=sum(Q(:,2,:))/dtz;
Qo(:,:)=sum(Q(:,3,:))/dtz;
Qp(:,:)=sum(Q(:,5,:)/dtz);
plot(T,Qz,T,Qd,T,Qo,T,Qp)

end

function plot_fild(x,y,z,c,Nl,X,Y,WXY,GY_subl,text,XYgy,a0,tn)


if Nl==1
    VORON3D
    hold on
    plot(WXY(a0==1,1),WXY(a0==1,2),'o','Color',[0 0 0],'MarkerFaceColor',[0.749019622802734 0 0.749019622802734]);
    plot(WXY(a0==-1,1),WXY(a0==-1,2),'o','Color',[0 0 0],'MarkerFaceColor',[0.249019622802734 0 0.949019622802734])
else
    figure1=figure(45);
    dx=max(x(:))-min(x(:));
    dy=max(y(:))-min(y(:));
    asx=dx/dy;

    axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[asx 1 1]);
    set(axes1,'PlotBoxAspectRatio',[asx 1 1])
    hold(axes1,'all');

    a=find(GY_subl(:,2)==-1);
    for l=1:size(a,1)
        if l<size(a,1)
            x_y(l)={GY_subl(a(l)+1:a(l+1)-1,:)};
        else
            x_y(l)={GY_subl(a(l)+1:end,:)};
        end
    end
    

    [~,ia,~] = intersect(WXY,[x(:,1),y(:,1)],'rows');
    WZ=z(ia,:);
    z_nan=WZ; z_nan(:)=nan;
        
    for i=1:Nl
        VORON3D(x(:,i),y(:,i),z(:,i),x_y{i},c(:,i))
        IN=inpolygon(WXY(:,1),WXY(:,2),x_y{i}(:,1),x_y{i}(:,2));
        z_nan(IN,i)=z(IN,i);
        hold on
    end;

%      wz=min(z(:))*ones(size(WXY,1),1);
%      wz(isnan(zii(:,end))==1)=z(isnan(zii(:,end))==1,end-1);
     plot3(repmat(WXY(:,1),1,2)',repmat(WXY(:,2),1,2)',[min(z_nan')',5+max(WZ(:))*ones(size(WXY,1),1)]','LineWidth',2,'Color',[0.749019622802734 0 0.749019622802734])
end;
 title(text);
 colorbar;
 hold off
end

function plot_crack_color(Nl,NT,SwC,CR_GRUP,XY,z,A2C)

if Nl==1
    for i1=1:size(NT,2)
        nt=NT{i1};
        SwC_nC=SwC(CR_GRUP==i1);
        unt=unique(nt);
        for i=1:size(nt,2)
            rr1=find(nt(1,i)==unt);
            rr2=find(nt(2,i)==unt);
            col_sw=(SwC_nC(rr1)+SwC_nC(rr2))/2;
            col_sw=((1-col_sw)*0.9+0.1);
            vx=[XY(nt(1,i),1),XY(nt(2,i),1)];
            vy=[XY(nt(1,i),2),XY(nt(2,i),2)];
            plot(vx,vy,'Color',[col_sw col_sw col_sw])
            
            hold on
        end;
    end;
else
    for i2=1:size(NT,2)
        nt=NT{i2};
        unt=unique(nt);
        sw=SwC(i2==CR_GRUP(:,1));
        a2c=A2C(:,i2==CR_GRUP(:,1));
        [r,c]=find(a2c);
    
        vxy=XY(r,:);
        vz=z(r);
        fz=zeros(size(z));
        fz(r)=1;
        [~,c1]=find(fz);
        na=size(sw,1);

        nl=max(c1)-min(c1)+1;
        vx=reshape(vxy(:,1),na/nl,nl);
        vy=reshape(vxy(:,2),na/nl,nl);
        vz1=reshape(vz,na/nl,nl);
        sw1=reshape(sw,na/nl,nl);
        surf(vx,vy,vz1,sw1,'LineStyle','none');
        %plot3(vx,vy,vz,'Color',[col_sw col_sw col_sw])
        hold on
    end;
end;
end

function VORON3D(x,y,z,BND,Sw)
xy=[x,y];
[V,C]=VoronoiLimit2(xy(:,1),xy(:,2),BND); 

Sw0=Sw;

% axes1 = axes('Parent',figure1,'YGrid','on','XGrid','on',...
%     'Position',[0.13 0.11 0.651609195402299 0.815],...
%     'FontSize',14);

for k=1:length(C)
    
    VertCell = V(C{k},:);
    ZZ=z(k)*ones(size(VertCell(:,1)));
    patch(VertCell(:,1),VertCell(:,2),ZZ,Sw0(k),'FaceColor','flat')
    
end
end