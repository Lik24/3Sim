CD=CD_rk_200_50;
c_off=0.98;
vp_off=4;

for i=1:size(CD,2)
 V0=CD{8,i};
 c(:,i)=CD{1,i};
 KIN(:,i)=CD{2,i}/V0;
 vp(:,i)=CD{3,i}/V0;
 
 c1(:,i)=c(:,i);
 c1(c(:,i)<=c_off,i)=c(c(:,i)<=c_off,i);
 c1(c(:,i)>c_off,i)=nan;
 r=find(isnan(c1(:,i))==0);
 
%  vp1(:,i)=vp(:,i);
%  vp1(vp(:,i)<=vp_off,i)=vp(vp(:,i)<=vp_off,i);
%  vp1(vp(:,i)>vp_off,i)=nan;
%  r=find(isnan(vp1(:,i))==0);
 
 kin1(:,i)=KIN(:,i);
 kin1(1:max(r),i)=KIN(1:max(r),i);
 kin1(max(r)+1:end,i)=nan;
 
 c1(:,i)=c(:,i);
 c1(1:max(r),i)=c(1:max(r),i);
 c1(max(r)+1:end,i)=nan;
 
 vp1(:,i)=vp(:,i);
 vp1(1:max(r),i)=vp(1:max(r),i);
 vp1(max(r)+1:end,i)=nan;
 
 c2(:,i)=c(:,i);
 c2(1:max(r),i)=c(1:max(r),i);
 c2(max(r)+1:end,i)=nan;
 T(i)=max(r);
 kn(i)=KIN(max(r),i);
 
 Vp(i)=vp(max(r),i);
end

cc2=c2;
for i=1:size(CD,2)
cc2(isnan(c2(:,i))==0,i)=smooth(c2(isnan(c2(:,i))==0,i));
end

plot_kin_obv(1:200, kin1(:,[1,2,4,6,7]), c2(:,[1,2,4,6,7]))
%plotyy((1:1800)/12,kin1(:,[1,3,5,7]),(1:1800)/12,c2(:,[1,3,5,7]))