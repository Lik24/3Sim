function Z_new=filterZ(Z,H,WXY)

nl=size(Z,2);

    x=repmat(WXY(:,1),nl,1);
    y=repmat(WXY(:,2),nl,1);
    x1=x;
    y1=y;
    k1=repmat([1:nl]',1,size(WXY,1));
    k1=reshape(k1',size(WXY,1)*nl,1);
    z=Z(:);
    k=k1;
    
    x(z==0)=[];
    y(z==0)=[];
    k(z==0)=[];    
    z(z==0)=[];
    
    
    F=scatteredInterpolant(x,y,k,z,'linear','nearest');
    V=F(x1,y1,k1);

    Z_new=reshape(V,size(WXY,1),nl);