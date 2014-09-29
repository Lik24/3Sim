function [gt,SS]=Tresh_Gor(n,gXY,NL)

A=MR_Prop(gXY,1);

for l=1:NL
    for i=1:n
        if i==1
            x=70;
            y=50;
            p=450;
        else
            x=700;
            y=500;
            p=300;
        end;
        XY=rand(10,2);
        XY(:,1)=x*XY(:,1)+p;
        XY(:,2)=y*XY(:,2)+p;
        a=convhull(XY(:,1),XY(:,2));
        XY=XY(a,:);
        S(:,i)=Area3(A,gXY,XY);
        gt(i,l)={find(S(:,i))};
    end;
        SS(l)={S};
end;
