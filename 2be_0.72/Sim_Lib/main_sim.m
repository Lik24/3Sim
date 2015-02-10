clear all
tic
as=[2,2,1];     %Степень фазовой w-o-g
aw=[1,1,10];    %Высота фазовой w-o-g
mu=[1,10,0.1];   %Вязкость w-o-g
Ns=2;          %Число скважин
Nl=1;           %Число слоёв
dt=5;           %Шаг по времени

[K,Mp,P,Sw,WXY,H,Z]=Sintetic(Ns,Nl);
[Pw,Uf]=Well_DATA(WXY,Z);

yxy(:,1)=[0:100:1000,1000*ones(1,9),1000:-100:0,0*ones(1,9)]';
yxy(:,2)=[0*ones(1,11),100:100:900,1000*ones(1,11),900:-100:100]';

%gXY=[rand(1000,2)*1000;yxy];
[gX,gY]=meshgrid(0:10:1000,0:10:1000);
gXY(:,1)=gX(:);
gXY(:,2)=gY(:);

[K,Mp,P,Sw,XY,H,Z,Won]=GridProp(K,Mp,P,Sw,WXY,H,Z,gXY);
K(:)=mean(K(:));

[A]=MR_Prop(XY,Nl);

[L,B,S,H1]=Geome3(A,XY,Z,H);
Wf=KWell(K,H,S,Won);

p=symrcm(A);
A=A(p,p);  L=L(p,p);
S=S(p,p);  B=B(p,p); 
H1=H1(p,p);
[r,c]=find(L);

XY=repmat(XY,Nl,1);
K=K(:); Mp=Mp(:); Sw=Sw(:); H=H(:); Z=Z(:); P=P(:); 
K=K(p); Mp=Mp(p); Sw=Sw(p); XY=XY(p,:); H=H(p); Z=Z(p); P=P(p); 

for i=1:size(Won,1)
   r1=find(Won(i)==p);
   Won(i)=r1;
end;

nc=size(XY,1);

figure(99),spy(A);
Pi(:,1)=P;
[Ke,He,dV]=KH2Mat(K,H,Mp,S,r,c); 

for t=1:2000
    [TL,TW]=Potok_MKT(L,B,Ke,Pi(:,t),Sw(:,t),H1,r,c,as,aw,mu);
    %[gTL,gTW]=Potok_MKT_GPU(L,B,Ke,Pi(:,t),Sw(:,t),H1,r,c,as,aw);
    [W1,W6]=Well_MKT(Wf,Won,Uf,Sw(:,t),aw,as,mu);
    A1=TL-sparse(Won,Won,W1,nc,nc);
    b1=sparse(Won,ones(1,size(Won,1)),-W1.*Pw',nc,1);
    Pi(:,t+1)=b1'/A1;
    A2=TW-sparse(Won,Won,W6,nc,nc);
    b2=sparse(Won,ones(1,size(Won,1)),W6.*Pw',nc,1);
    Sw(:,t+1)=Sw(:,t)+dt*(A2*Pi(:,t+1)+b2)./dV;
    Sw(:,t+1)=Sw(:,t+1).*(Sw(:,t+1)>=0).*(Sw(:,t+1)<=1)+(Sw(:,t+1)>1);
end;

Q=QBild(W1,W6,b1,b2,Pi,Uf,Won,nc);

VZL(XY,K,WXY,Z,Pi,Sw,Nl,p);

toc