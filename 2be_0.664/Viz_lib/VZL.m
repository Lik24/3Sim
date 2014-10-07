function VZL(XY,K,WXY,WZ,P,Sw,T,Cp,Nl,pp,Q,PXY,gt)

 %cmap1=summer;
 %colormap(cmap1);


P=P(1:size(XY,1),:);
T=T(1:size(XY,1),:);
%K(pp)=K;
XY(pp,1)=XY(:,1);
XY(pp,2)=XY(:,2);
P(pp,:)=P;
Sw(pp,:)=Sw;
Cp(pp,:)=Cp;
T(pp,:)=T;
nxy=size(WZ,1)/Nl;

for l=1:Nl
    
Z(nxy*(l-1)+1:nxy*l)=(-l+1)*20;
end;

x=reshape(XY(:,1),size(WZ,1)/Nl,Nl);
y=reshape(XY(:,2),size(WZ,1)/Nl,Nl);
z=reshape(Z,size(WZ,1)/Nl,Nl);
log10k=reshape((K/8.64),size(WZ,1)/Nl,Nl);
p=reshape(P(:,end),size(WZ,1)/Nl,Nl);
sw=reshape(Sw(:,end),size(WZ,1)/Nl,Nl);
cp=reshape(Cp(:,end),size(WZ,1)/Nl,Nl);
tt=reshape(T(:,end),size(WZ,1)/Nl,Nl);

mx(1)=min(XY(:,1));
mx(2)=max(XY(:,1));
my(1)=min(XY(:,2));
my(2)=max(XY(:,2));
[X,Y]=meshgrid(mx(1):5:mx(2),my(1):5:my(2));
 
figure(98),subplot(2,4,1);
plot_fild(x,y,z,p,Nl,X,Y,WXY,'Пластовое давление') % 
hold on

if Nl==1
    pxy=PXY{1};
    for i1=1:size(pxy,2)
        pX=pxy{1,i1};
        pY=pxy{2,i1};
        
        for i=1:size(pX,1)
            plot(pX(i,:),pY(i,:),'Color',[0 0 i1/size(pxy,2)])
            hold on
        end;
    end;
else
    for i2=1:Nl
      pXY=PXY{i2};
       for j1=1:size(pXY,1)
            pxy=pXY{j1};
           for i1=1:size(pxy,2)
            pX=pxy{1,i1};
            pY=pxy{2,i1};
            pZ=mean(z(:,i2))*ones(size(pY));
            if isempty(pX)==0
                plot3(pX',pY',pZ','Color',[0 0 i2/size(PXY,2)])
            end;
           end;
            hold on
        end;
    end;
end;
hold off

figure(98),subplot(2,4,2);
plot_fild(x,y,z,sw,Nl,X,Y,WXY,'Водонасыщенность') % 
hold on

if Nl==1
    pxy=PXY{1};
    for i1=1:size(pxy,2)
        pX=pxy{1,i1};
        pY=pxy{2,i1};
        
        for i=1:size(pX,1)
            plot(pX(i,:),pY(i,:),'Color',[0 0 i1/size(pxy,2)])
            hold on
        end;
    end;
else
    for i2=1:Nl
      pXY=PXY{i2};
       for j1=1:size(pXY,1)
            pxy=pXY{j1};
           for i1=1:size(pxy,2)
            pX=pxy{1,i1};
            pY=pxy{2,i1};
            pZ=mean(z(:,i2))*ones(size(pY));
            if isempty(pX)==0
                plot3(pX',pY',pZ','Color',[0 0 i2/size(PXY,2)])
            end;
           end;
            hold on
        end;
    end;
end;
hold off


figure(98),subplot(2,4,4);
plot_fild(x,y,z,cp,Nl,X,Y,WXY,'Концентрация') % 

figure(98),subplot(2,4,3);
plot_fild(x,y,z,tt,Nl,X,Y,WXY,'Температура') % 

figure(98),subplot(2,4,5);
plot_fild(x,y,z,log10k,Nl,X,Y,WXY,'Проницаемость') % 



% 
% 
% for i2=1:size(gt,2)
%    for i1=1:size(gt,1)
%     a=gt{i1,i2};   
%     if isempty(a)==0
%     b=convhull(XY(a,1),XY(a,2));
%     plot(XY(a(b),1),XY(a(b),2),'m--')
%     end;
%    end
% end;
% hold off
% 
% subplot(2,4,2);
% contourf(X,Y,mean(C,3))
% 


 subplot(2,4,8);
 T=1:size(Q,3);
% 
Qz(:,:)=sum(Q(:,1,:));
Qd(:,:)=sum(Q(:,2,:));
Qo(:,:)=sum(Q(:,3,:));
Qp(:,:)=sum(Q(:,5,:));
plot(T,Qz,T,Qd,T,Qo,T,Qp)
end

function plot_fild(x,y,z,c,Nl,X,Y,WXY,text)
if Nl==1
    F=scatteredInterpolant(x(:,1),y(:,1),c(:,1),'nearest','linear');
    K=F(X,Y);
    contourf(mean(X,3),mean(Y,3),mean(K,3),'LineStyle','none');
    hold on
    plot(WXY(:,1),WXY(:,2),'o','Color',[0 0 0],'MarkerFaceColor',[0.749019622802734 0 0.749019622802734])

else
    %wz=reshape(WZ,size(WZ,1)/Nl,Nl);
    for i=1:Nl
        F=scatteredInterpolant(x(:,i),y(:,i),c(:,i),'nearest','none');
        K(:,:,i)=F(X,Y);
        F=scatteredInterpolant(x(:,i),y(:,i),z(:,i),'nearest','none');
        Z1(:,:,i)=F(X,Y);
    end;
    
    for i=1:Nl
        surf(X,Y,Z1(:,:,i),K(:,:,i),'LineStyle','none')
        hold on
    end;

     plot3(repmat(WXY(:,1),1,2)',repmat(WXY(:,2),1,2)',[min(z(:))*ones(size(WXY,1),1),5+max(z(:))*ones(size(WXY,1),1)]','LineWidth',2,'Color',[0.749019622802734 0 0.749019622802734])
end;
 title(text);
 colorbar;
 hold off
end
