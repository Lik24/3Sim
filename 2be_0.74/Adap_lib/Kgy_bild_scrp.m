load('E:\Облако\Юни-Конкорд\2014\Бузачи_3Sim\2be_0.664\Adap_lib\ADP_ReZ2.mat','GY_Kxy')
nw=sum(u1(u1==1));
r=find(u1>0);

BndXY=DATA.BndXY;
K=zeros(177,nw);
K(BndXY==1,:)=Kxy(:,r);


for i=1:size(K,1)
  K_new(i,1)=mean([GY_Kxy(i,end),K(i,K(i,:)~=0)]);
  if sum(K(i,:)~=0)~=0
      K_new(i,1)=mean(K(i,K(i,:)~=0));
  end;
end;