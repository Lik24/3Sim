function VZL_pik(DATA,WXY,P,Sw,Sw0,Nl,pp,Q,SwC,CR_GRUP,pc,NT,XYgy,a0,i)
j=i*3;
XY=DATA.XY;
WZ=DATA.gZ;
ka=DATA.ka;
K=DATA.gKX;

XY=repmat(XY,Nl,1);

P1=P(1:size(pp,2),:);


P1(pp,:)=P1;
Sw1(pp,:)=Sw;
SwC(pc,:)=SwC;
KIN1=1-(1-Sw1(:,end))./(1-Sw0);

P=zeros(size(ka,1),size(P,2));
Sw=zeros(size(ka,1),size(P,2));
KIN=zeros(size(ka,1),1);

P(ka==1,:)=P1;
Sw(ka==1,:)=Sw1;
KIN(ka==1,:)=KIN1;

P(ka==0,:)=nan;
Sw(ka==0,:)=nan;
KIN(ka==0,:)=nan;
nxy=size(WZ,1)/Nl;
%[A]=MR_Prop(XY,1);
%A=A.*(A>0);
for l=1:Nl
Z(nxy*(l-1)+1:nxy*l)=(-l+1)*20;
end;

x=reshape(XY(:,1),size(WZ,1)/Nl,Nl);
y=reshape(XY(:,2),size(WZ,1)/Nl,Nl);
z=reshape(Z,size(WZ,1)/Nl,Nl);

p=reshape(P(:,end),size(WZ,1)/Nl,Nl);
sw=reshape(Sw(:,end),size(WZ,1)/Nl,Nl);
kin=reshape(KIN(:,end),size(WZ,1)/Nl,Nl);

mx(1)=min(XY(:,1));
mx(2)=max(XY(:,1));
my(1)=min(XY(:,2));
my(2)=max(XY(:,2));
[X,Y]=meshgrid(mx(1):5:mx(2),my(1):5:my(2));
 
figure1=figure(j-2);
plot_fild(x,y,z,p,Nl,X,Y,WXY,'Пластовое давление',XYgy,a0,j-2,i) % 
saveas(figure1,strcat('pic_Вар№',num2str(i),'_дав.jpg'))
% hold on
% plot_crack_color(Nl,NT,SwC,CR_GRUP,XY,z);
% hold off

figure2=figure(j-1);
plot_fild(x,y,z,1-sw,Nl,X,Y,WXY,'Нефтенасыщенность',XYgy,a0,j-1,i) % 
saveas(figure2,strcat('pic_Вар№',num2str(i),'_нас.jpg'))
% hold on
% plot_crack_color(Nl,NT,SwC,CR_GRUP,XY,z);
% hold off
figure3=figure(j);
plot_fild(x,y,z,kin,Nl,X,Y,WXY,'КИН',XYgy,a0,j,i) % 
saveas(figure3,strcat('pic_Вар№',num2str(i),'_кин.jpg'))

%  subplot(2,4,8);
%  T=1:size(Q,3);
% % 
% Qz(:,:)=sum(Q(:,1,:));
% Qd(:,:)=sum(Q(:,2,:));
% Qo(:,:)=sum(Q(:,3,:));
% Qp(:,:)=sum(Q(:,5,:));
% plot(T,Qz,T,Qd,T,Qo,T,Qp)

end

function plot_fild(x,y,z,c,Nl,X,Y,WXY,text,XYgy,a0,jj,j1)
dx=max(X(1,:))-min(X(1,:));
dy=max(Y(:,1))-min(Y(:,1));
asx=dx/dy;

axes1 = axes('Parent',figure(jj),'PlotBoxAspectRatio',[asx 1 1]);
set(axes1,'PlotBoxAspectRatio',[asx 1 1])
hold(axes1,'all');

if Nl==1
    F=scatteredInterpolant(x(:,1),y(:,1),c(:,1),'nearest','none');
    K=F(X,Y);
    [IN,ON]=inpolygon(X(:),Y(:),XYgy(:,1),XYgy(:,2));
    K([IN+ON]==0)=NaN;
    contourf(mean(X,3),mean(Y,3),mean(K,3),'LineStyle','none');
    hold on
    plot(WXY(a0==1,1),WXY(a0==1,2),'o','Color',[0 0 0],'MarkerFaceColor',[0.749019622802734 0 0.749019622802734]);
    plot(WXY(a0==-1,1),WXY(a0==-1,2),'o','Color',[0 0 0],'MarkerFaceColor',[0.249019622802734 0 0.949019622802734])
    plot(WXY(a0==0,1),WXY(a0==0,2),'o','Color',[0 0 0],'MarkerFaceColor',[1 1 1])
    grid on
else
    %wz=reshape(WZ,size(WZ,1)/Nl,Nl);
    for i=1:Nl
        F=scatteredInterpolant(x(:,i),y(:,i),c(:,i),'nearest','none');
        K(:,:,i)=F(X,Y);
        [IN,ON]=inpolygon(X(:),Y(:),XYgy(:,1),XYgy(:,2));
        KK=K(:,:,i);
        KK([IN+ON]==0)=NaN;
        K(:,:,i)=KK;
        F=scatteredInterpolant(x(:,i),y(:,i),z(:,i),'nearest','none');
        Z1(:,:,i)=F(X,Y);
    end;
    
    for i=1:Nl
        surf(X,Y,Z1(:,:,i),K(:,:,i),'LineStyle','none')
        hold on
    end;

     plot3(repmat(WXY(:,1),1,2)',repmat(WXY(:,2),1,2)',[min(z(:))*ones(size(WXY,1),1),5+max(z(:))*ones(size(WXY,1),1)]','LineWidth',2,'Color',[0.749019622802734 0 0.749019622802734])
end;
 title(strcat('Вариант №',num2str(j1)));
 legend(text,'Доб. скв.','Наг. скв.','Откл. скв.')
 colorbar;
 hold off
end

function plot_crack_color(Nl,NT,SwC,CR_GRUP,XY,z)

if Nl==1
    Nt=NT{1};
    SwC_L=SwC(CR_GRUP(:,2)==Nl,end);
    cr_grup=CR_GRUP(CR_GRUP(:,2)==Nl,1);
    
    for i1=1:size(Nt,2)
        nt=Nt{i1};
        SwC_nC=SwC_L(cr_grup==i1);
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
    for i2=1:Nl
        Nt=NT{i2};
        SwC_L=SwC(CR_GRUP(:,2)==i2,end);
        cr_grup=CR_GRUP(CR_GRUP(:,2)==i2,1);
        
        for i1=1:size(Nt,2)
            nt=Nt{i1};
            SwC_nC=SwC_L(cr_grup==i1);
            unt=unique(nt);
            for i=1:size(nt,2)
                rr1=find(nt(1,i)==unt);
                rr2=find(nt(2,i)==unt);
                col_sw=(SwC_nC(rr1)+SwC_nC(rr2))/2;
                col_sw=((1-col_sw)*0.9+0.1);
                vx=[XY(nt(1,i),1),XY(nt(2,i),1)];
                vy=[XY(nt(1,i),2),XY(nt(2,i),2)];
                vz=[z(nt(1,i),i2),z(nt(2,i),i2)];
                
                plot3(vx,vy,vz,'Color',[col_sw col_sw col_sw])
                
                hold on
            end;
        end;
    end;
end;
end