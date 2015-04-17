function Wf=KWell_Horiz(Wf,KX,KY,KZ,H,S,L,B,Won,r,c,Perf,SDoly,r0,XY,Nw,Nl)
Won=Won(1:size(Won,1)/Nl,:);

if sum(Won(:,2)==1)>0
    WGOR=Won(Won(:,2)==1,3);
    uWN=unique(WGOR);
    for i=uWN'
        temp_w=Won(i==Won(:,3),:);
        prd=temp_w(:,4);
        [B,I]=sort(prd);
        temp_w=temp_w(I,:);
        l=zeros(size(temp_w,1),1);
        for j=1:size(temp_w,1)
            if j==1
                l(j)=L(temp_w(j,1),temp_w(j+1,1));
            elseif j==size(temp_w,1)
                l(j)=L(temp_w(j-1,1),temp_w(j,1));
            else
                l(j)=L(temp_w(j-1,1),temp_w(j,1))+L(temp_w(j,1),temp_w(j+1,1));
            end
        end
        kx=KX(temp_w(:,1),:);
        kz=KZ(temp_w(:,1),:);
        dx2=repmat(l.^2,1,Nl);
        dz2=H(temp_w(:,1),:).^2;
        Rh=0.28*(dx2.*(kz./kx).^0.5+dz2.*(kx./kz).^0.5).^0.5./((kz./kx).^0.25+(kx./kz).^0.25);
        Rh(I,:)=Rh;
        Rhh(i==Won(:,3),:)=Rh;
        L1=l;
        L1(I,:)=l;
        LL(i==Won(:,3),:)=L1;
    end
    
    v1=Won(:,2)==1;
    Wcof=zeros(size(v1,1),Nl);
    for l=1:Nl
        Wcof(v1,l)=(KX(Won(v1,1),l).*KZ(Won(v1,1),l)).^0.5.*LL(v1)./log(Rhh(v1)/r0).*Perf(v1,l).*SDoly(v1);
    end;
    Wf(repmat(v1,Nl,1))=2*pi*Wcof(repmat(v1,Nl,1));
end