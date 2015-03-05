function [ng]=Tresh_Gor(on,XY,Nl)
ng=cell(Nl,1);
if on==1
    n=size(XY,1);
    for l=1:Nl
        LN=n*(l-1);
        for i=1:1
            xy=XY(1:10,:);
            a=convhull(xy(:,1),xy(:,2));
            nd_sl={a+LN};
        end
        ng(l)={nd_sl};
    end
else
     ng0={[]};
     ng(:)={ng0};
end