T=1:150;
for i=[2,4,6 ]%size(CD1,2)
 Q=CD1{1,i};
 c=1-Q(:,1)./Q(:,2);
 C(i)=c(end);
 kin=Q(:,4)/v0(i);
 KIN(i)=kin(end);
 kin1=kin(c<0.95);
 KIN(i)=kin1(end);
 plotyy(T,Q(:,4)/v0(i),T,1-Q(:,1)./Q(:,2))
 hold on
end