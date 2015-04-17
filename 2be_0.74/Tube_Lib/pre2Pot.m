function DZA=pre2Pot(dZ,n,RC,rc_in_h,rc_in_hd)

for i=1:4
    if i==1
       r=rc_in_h(:,1);
       c=rc_in_h(:,2);
    elseif i==2
       r=RC.Cr2;
       c=RC.Cc2;
    elseif i==3
       r=RC.Gr2;
       c=RC.Gc2;
    else
       r=rc_in_hd(:,1);
       c=rc_in_hd(:,2);
    end
    if n(i)~=0    
    DZA(i).A=muf(dZ(i,:),r,c,n(i));
    else
    DZA(i).A=zeros(0,2); 
    end
end

end

function A=muf(DZ,r,c,n)
    dZW=DZ{1};  dZO=DZ{2};  dZG=DZ{3};
    v=r+(c-1)*n;
    A(:,1)=dZW(v);
    A(:,2)=dZO(v);
end