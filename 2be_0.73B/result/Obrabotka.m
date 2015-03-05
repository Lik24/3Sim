for i=1:5
kin1(:,i)=CD5{i}(:,4);
kin2(:,i)=CD6{i}(:,4);
kin3(:,i)=CD7{i}(:,4);

c1(:,i)=CD5{i}(:,1);
c2(:,i)=CD6{i}(:,1);
c3(:,i)=CD7{i}(:,1);

ki1=kin1(c1(:,i)<=0.98,i);
ki2=kin2(c2(:,i)<=0.98,i);
ki3=kin3(c3(:,i)<=0.98,i);

k1(i)=max(ki1);
k2(i)=max(ki2);
k3(i)=max(ki3);

t1(i)=find(ki1==k1(i));
t2(i)=find(ki2==k2(i));
t3(i)=find(ki3==k3(i));
end

