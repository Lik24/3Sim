

 dir1=XY1(2,:)-XY1(1,:);
 dir2=XY2(2,:)-XY2(1,:);
 
 a1=-dir1(2);
 b1=dir1(1);
 d1=-(a1*XY1(1,1)+b1*XY1(1,2));
 
 a2=-dir2(2);
 b2=dir2(1);
 d2=-(a2*XY2(1,1)+b2*XY2(1,2));
 
 seg1_l2_s=a2*XY1(1,1)+b2*XY1(1,2)+d2;
 seg1_l2_e=a2*XY1(2,1)+b2*XY1(2,2)+d2;
 seg2_l1_s=a1*XY2(1,1)+b2*XY2(1,2)+d1;
 seg2_l1_e=a1*XY2(2,1)+b2*XY2(2,2)+d1;
 
 fl=1;
 if seg1_l2_s*seg1_l2_e>=0 || seg2_l1_s*seg2_l1_e>=0
    fl=0;
 end
 u=seg1_l2_s/(seg1_l2_s-seg1_l2_e);
 A=XY1(1,:)+u*dir1;
