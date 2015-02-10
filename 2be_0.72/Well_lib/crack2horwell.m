function [WXY,Won,g_cr,de]=crack2horwell(g_cr,WXY,Won,drob)
de=[];
if sum(Won(:,2))~=0
    won=Won(Won(:,2)==1,:);
    wxy=WXY(Won(:,2)==1,:);
    new_XY=zeros(0,2);
    new_Won=zeros(0,3);
    for i=1:size(g_cr,1)
        XY1=g_cr{i}(:,1:2);
        Z=g_cr{i}(:,3);
        WN=unique(won(:,1));
        new_XY1=zeros(0,2);
        for k=1:size(WN,1)
            xy=wxy(WN(k)==won(:,1),:);
            por=won(WN(k)==won(:,1),3);
            [B,I]=sort(por);
            por=por(I);
            xy=xy(I,:);
            A=[];
            fl=[];
            for i1=1:size(xy,1)-1
                XY2=xy(i1:i1+1,:);
                [A(i1,:),fl(i1),flm(i1),flp(i1)]=find_poin(XY1,XY2);
            end
            new_xy=A(fl==1,:);
            if sum(fl)>0
                new_won=[WN(k)*ones(sum(fl),1),ones(sum(fl),1),por(fl==1)+0.5.*(flp(fl==1)==1)'-0.5.*(flm(fl==1)==1)'];
            else
                new_won=zeros(0,3);
            end;
            [new_xy,ia,ic]=unique(new_xy,'rows');
            new_won=new_won(ia,:);
            new_XY=[new_XY;new_xy];
            new_XY1=[new_XY1;new_xy];
            new_Won=[new_Won;new_won];
        end
        g_cr{i}=[XY1(1,:),Z(1);new_XY1,mean(Z);XY1(2,:),Z(2)];
    end;
    
    for i1=1:size(new_XY,1)
      R=((WXY(:,1)-new_XY(i1,1)).^2+(WXY(:,2)-new_XY(i1,2)).^2).^0.5;
      de=[de;find(R<drob/2)];
    end
      WXY(de,:)=[];
      Won(de,:)=[];
      
    WXY=[WXY;new_XY];
    Won=[Won;new_Won];
end

end

function [A,fl,flm,flp]=find_poin(XY1,XY2)
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
 seg2_l1_s=a1*XY2(1,1)+b1*XY2(1,2)+d1;
 seg2_l1_e=a1*XY2(2,1)+b1*XY2(2,2)+d1;
 
 fl=1;
 flm=seg2_l1_s==0;
 flp=seg2_l1_e==0;
    
 if seg1_l2_s*seg1_l2_e>0 || seg2_l1_s*seg2_l1_e>0
    fl=0;
 end
 u=seg1_l2_s/(seg1_l2_s-seg1_l2_e);
 A=XY1(1,:)+u*dir1;
end