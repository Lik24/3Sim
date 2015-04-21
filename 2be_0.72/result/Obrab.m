
% for i=[3]%1:size(CD,2)
%     for j=1:5
%         kin=CD{j,i}(:,4);
%         %kin=kin/(500*500*40*0.3*0.075);
%         c=CD{j,i}(:,1);
%         vp=CD{j,i}(:,3);
%         ki=kin(c<0.98);
%         vpi=vp(c<0.98);
%         KIN(i,j)=max(ki);
%         Vp(i,j)=max(vpi)/(500*500*40*0.3);
%         r=find(kin==max(ki));
%         T(i,j)=r;
%     end;
% end


for i=[2,3,6]%1:size(CD,2)
    for j=1:4
        kin=CD{j,i}(:,4);
        %kin=kin/(500*500*40*0.3*0.075);
        c=CD{j,i}(:,1);
        vp=CD{j,i}(:,3);
        ki=kin(c<0.98);
        vpi=vp(c<0.98);
        KIN(i,j)=max(ki);
        Vp(i,j)=max(vpi)/(500*500*40*0.3);
        r=find(kin==max(ki));
        T(i,j)=r;
    end;
end