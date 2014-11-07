
for i=1:size(SD,2)
CD=SD{i};
Q=CD{1};
Pi=CD{2};
Sw=CD{3};
p=CD{4};
PpW=CD{5};
dV0=CD{6};
DATA=CD{7};
WXY=CD{8};
WData=CD{9};
Sw0=CD{10};
PR=CD{11};
XYc=CD{12};
Pw=CD{13};

Qo(:,i)=sum(Q(:,3,:));
Ql(:,i)=sum(Q(:,2,:));
Qz(:,i)=sum(Q(:,1,:));

QL(:,:)=Q(:,2,:);
QZ(:,:)=Q(:,1,:);

C(:,i)=1-Qo(:,i)./Ql(:,i);
V0=sum(dV0.*(1-Sw0));
% [(C(:,i)<0.98),(isnan(C(:,i))==0),(C(:,i)<0.98).*(isnan(C(:,i)==0))]'
uf=WData.Uf(:,(C(:,i)<0.98).*(isnan(C(:,i))==0)==1);
VZL_pik(DATA,WXY,Pi,Sw,Sw0,PR.Nl,p,Q,[],[],[],[],XYc,uf(:,end-10),i);

Nw_dob(:,i)=sum(WData.Uf==1,1);
Nw_zak(:,i)=sum(WData.Uf==-1,1);
Ppw(:,i)=mean(PpW,1);
W1=(QL+QZ)./(PpW-Pw);
 for t=1:size(Pw,2) 
    PW_d(t,i)=mean(Pw(WData.Uf(:,t)==1,t),1);
    PW_i(t,i)=mean(Pw(WData.Uf(:,t)==-1,t),1);
    W11=W1(WData.Uf(:,t)==1,t);
    W1_d(t,i)=mean(W11(isnan(W11)==0),1);
    W12=W1(WData.Uf(:,t)==-1,t);
    W1_i(t,i)=mean(W12(isnan(W12)==0),1);
 end

end;

sQo=cumsum(Qo,1);
sQl=cumsum(Ql,1);
sQz=-cumsum(Qz,1);
Vp=sQz/sum(dV0);

KIN=sQo/V0;
T=1:size(KIN,1);
T=T/365;
KIN(C>0.98)=nan;%
figure(i*3+1)
plot(T,KIN);
legend('������� �1','������� �2','������� �3');
grid on
T=repmat(T,i,1);
lT=T(C'<=0.98);
mT=max(lT(:));
xlim([0 ceil(mT)]);
xlabel('������, ���.');
ylabel('���, �.��.');
saveas(figure(i*3+1),strcat('pic_���.jpg'))

sQo(C>0.98)=nan;%
figure(i*3+2)
plot(T',sQo);
legend('������� �1','������� �2','������� �3');
grid on
lT=T(C'<=0.98);
mT=max(lT(:));
xlim([0 ceil(mT)]);
xlabel('������, ���.');
ylabel('����������� �����, �^3');
saveas(figure(i*3+2),strcat('pic_�����. �����.jpg'))

figure(i*3+3)
plot(Vp,KIN);
legend('������� �1','������� �2','������� �3');
grid on
lVp=Vp(C<=0.98);
mVp=max(lVp(:));
xlim([0 ceil(mVp)]);
xlabel('����������� ������� �����, ���. ��.');
ylabel('���, �.��.');
saveas(figure(i*3+3),strcat('pic_���_VP.jpg'))

sQz(C>0.98)=nan;%
Qz(C>0.98)=nan;%
figure(i*3+4)
plot(T',sQz./sQl);
hold on
plot(T',-Qz./Ql,'--');
legend('������� �1 �����.','������� �2 �����.','������� �3 �����.','������� �1 ���.','������� �2 ���.','������� �3 ���.');
grid on
lT=T(C'<=0.98);
mT=max(lT(:));
xlim([0 ceil(mT)]);
xlabel('������, ���.');
ylabel('�����������, �.��.');
saveas(figure(i*3+4),strcat('pic_����.jpg'))

C(C>0.98)=nan;%
figure(i*3+5)
plot(T',C);
legend('������� �1','������� �2','������� �3');
grid on
lT=T(C'<=0.98);
mT=max(lT(:));
xlim([0 ceil(mT)]);
xlabel('������, ���.');
ylabel('������������, �.��.');
saveas(figure(i*3+5),strcat('pic_���.jpg'))