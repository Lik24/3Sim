function Wf=KWell_3(K,H,S,L,B,Won,r,c,Doly,r0,Con,WXY)
Perf=1;
Con=[Con;[WXY,(1:size(WXY,1))']];
C1=Con(:,1); C2=Con(:,2); C3=Con(:,3);

for i=1:size(Won,1);
    x=r(c==Won(i));
    y=c(c==Won(i));
    Con1=[C1(Won(i)==Con(:,3)),C2(Won(i)==Con(:,3)),C3(Won(i)==Con(:,3))];
    C=convhull([Con1(:,1),Con1(:,2)]);
    tet(i,1)=sumygl([Con1(C,1),Con1(C,2)],WXY(i,1),WXY(i,2));
    for j=1:size(x,1)
        l1(j)=L(x(j),y(j));
        b1(j)=B(x(j),y(j));
    end;
        Rk(i,1)=exp((sum(b1./l1.*log(l1),2)-tet(i,1))/sum(b1./l1,2));
        clear b1 l1 Con1
end;

Wcof=K(Won).*H(Won)./log(Rk/r0);
Wf=2*pi*Wcof;
end

function T=sumygl(A,x1,y1)
T=0;
for i=1:size(A,1)-1
    if ((x1~=A(i,1))+(y1~=A(i,2))).*((x1~=A(i+1,1))+(y1~=A(i+1,2)))
        T=T+ygol(x1,y1,A(i,1),A(i,2),A(i+1,1),A(i+1,2));
    end;
end;
end

function TET=ygol(x1,y1,x2,y2,x3,y3)
xa=x2-x1;
ya=y2-y1;
xb=x3-x1;
yb=y3-y1;
TET=acos((xa*xb+ya*yb)/(((xa^2+ya^2)*(xb^2+yb^2))^0.5));
TET=real(TET);
end