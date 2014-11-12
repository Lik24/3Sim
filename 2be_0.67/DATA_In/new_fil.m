nl=size(Z,2);

mx(1)=min(WXY(:,1));
mx(2)=max(WXY(:,1));
my(1)=min(WXY(:,2));
my(2)=max(WXY(:,2));
 mz(1)=min(Z(:));
 mz(2)=max(Z(:));
[X,Y,Z1]=meshgrid(mx(1):50:mx(2),my(1):50:my(2),mz(1):5:mz(2));

%for l=1:nl
    x=repmat(WXY(:,1),nl,1);
    y=repmat(WXY(:,2),nl,1);
   % x=WXY(:,1);
   % y=WXY(:,2);
    x1=x;
    y1=y;
    k1=repmat([1:nl]',1,size(WXY,1));
    k1=reshape(k1',size(WXY,1)*nl,1);
    z=Z(:);
    %k=k1;
    k1=Z(:,1);
    k1(k1==0)=mean(k1(k1~=0));
    
    k2=Z(:,2);
    k2(k2==0)=mean(k2(k2~=0));
    
    k3=Z(:,3);
    k3(k3==0)=mean(k3(k3~=0));
    k=[k1;k2;k3];
    %x(z==0)=[];
    %y(z==0)=[];
    %k(z==0)=[];    
    z(z==0)=nan;
    
    F=scatteredInterpolant(x,y,k,z,'linear','nearest');
    Vz=F(X,Y,Z1);
    yw = inpaintn(Vz);
    
 for l=1:nl
    figure(l),subplot(1,2,1),contourf(X(:,:,l),Y(:,:,l),Vz(:,:,l))
    figure(l),subplot(1,2,2),contourf(X(:,:,l),Y(:,:,l),yw(:,:,l))    
 end
    Vz1=yw;
%end

Z3=Z;
f0=find(sum(H,2)==0);
[A]=MR_Prop(WXY,1);


for l=1:nl
    vz=Vz1(:,:,l);
    F=scatteredInterpolant(X(:),Y(:),vz(:),'linear','nearest');
    for i=1:size(f0,1)
        Z3(f0(i),l)=F(WXY(f0(i),1),WXY(f0(i),2));
    end
end
%  
%  K1=K;
%  Mp1=Mp;
%  Sw1=Sw;
%  H1=H;
%  Hk1=Hk;
%  Z1=Z;
%  
% f0=find(sum(H,2)==0); 
% [A]=MR_Prop(WXY,1);
% 
% for i=1:size(f0,1)
%     A1=A(:,f0(i));
%     A1(f0)=0;
%     c=find(A1);
%     if isempty(c)==0
%         for l=1:3
%             K1(f0(i),l)=mean(K(c,l));
%             Mp1(f0(i),l)=mean(Mp(c,l));
%             Sw1(f0(i),l)=mean(Sw(c,l));
%             H1(f0(i),l)=mean(H(c,l));
%             Hk1(f0(i),l)=mean(Hk(c,l));
%             Z1(f0(i),l)=mean(Z(c,l));
%         end;
%     else
%         f1=[f1,f0(i)];
%     end
% end;
