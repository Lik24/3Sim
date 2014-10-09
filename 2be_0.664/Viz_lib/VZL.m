function VZL(XY,K,WXY,WZ,P,Sw,T,Cp,Nl,pp,Q,SwC,CR_GRUP,pc,NT)

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
SwC(pc,:)=SwC;

nxy=size(WZ,1)/Nl;
%[A]=MR_Prop(XY,1);
%A=A.*(A>0);
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
plot_crack_color(Nl,NT,SwC,CR_GRUP,XY,z);
hold off

figure(98),subplot(2,4,2);
plot_fild(x,y,z,sw,Nl,X,Y,WXY,'Водонасыщенность') % 
hold on
plot_crack_color(Nl,NT,SwC,CR_GRUP,XY,z);
hold off


figure(98),subplot(2,4,4);
plot_fild(x,y,z,cp,Nl,X,Y,WXY,'Концентрация') % 

figure(98),subplot(2,4,3);
plot_fild(x,y,z,tt,Nl,X,Y,WXY,'Температура') % 

figure(98),subplot(2,4,5);
plot_fild(x,y,z,log10k,Nl,X,Y,WXY,'Проницаемость') % 



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

function plot_crack_color(Nl,NT,SwC,CR_GRUP,XY,z)
size(z)
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
                vz=[z(nt(1,i),2),z(nt(2,i),2)];
                
                plot3(vx,vy,vz,'Color',[col_sw col_sw col_sw])
                
                hold on
            end;
        end;
    end;
end;
end