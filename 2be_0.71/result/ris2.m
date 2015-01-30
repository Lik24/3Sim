%load('matlab16.12.mat');

for i=1:size(SD,2)
  CD=SD{i};
  Q=CD{1};
  Qo(:,i)=sum(Q(:,3,:));
  Ql(:,i)=sum(Q(:,2,:));
end;

for i=1:size(Qo,1)
 sQo(i,:)=sum(Qo(1:i,:),1);
 sQl(i,:)=sum(Ql(1:i,:),1);
end;

V=250*250*10*0.2;
KIN=sQo/V;
c=1-Qo./Ql;
T=1:365*50;
T=T/365;
Per=[10,5,2,1,0.5,0.1];
v=[1,2:8];
mS=18;

figure(1),axes('Parent',figure(1),'FontSize',mS);
plot(T,c(:,v),'LineWidth',2);
xlabel('Период, год.','FontSize',mS);    ylabel('Обводнённость, д.ед.','FontSize',mS);
grid on

figure(2),axes('Parent',figure(2),'FontSize',mS);
plot(T,KIN(:,v),'LineWidth',2);
xlabel('Период, год.','FontSize',mS);    ylabel('КИН, д.ед.','FontSize',mS);
grid on

figure(3),axes('Parent',figure(3),'FontSize',mS);
plot(T,Qo(:,v),'LineWidth',2);
xlabel('Период, год.','FontSize',mS);    ylabel('Дебит нефти, м^3.','FontSize',mS);
grid on

figure(4),axes('Parent',figure(4),'FontSize',mS);
plot(sQl(:,v)/V,KIN(:,v),'LineWidth',2);
xlabel('Прокачанный поровый объём, д. ед.','FontSize',mS);    ylabel('КИН, д.ед.','FontSize',mS);
grid on

figure(5),axes('Parent',figure(5),'FontSize',mS);
bar(c(end,v));
xlabel('Период, год.','FontSize',mS);    ylabel('Обводнённость, д.ед.','FontSize',mS);
grid on

figure(6),axes('Parent',figure(6),'FontSize',14);
bar(KIN(end,v));
xlabel('Период, год.','FontSize',mS);    ylabel('КИН, д.ед.','FontSize',mS);
grid on