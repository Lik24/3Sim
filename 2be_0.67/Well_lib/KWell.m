function Wf=KWell(K,H,S,L,B,Won_all,r,c,Perf,SDoly,r0,XY,Nw,Nl)

Won=Won_all(1:Nw);
WXY=XY(Won,:);

tet=2*pi*ones(size(WXY,1),1);
V=convhull(WXY);

for i=1:size(Won,1);

    x=r(c==Won(i));
    y=c(c==Won(i));

    if ~isempty(find(V==i))
        Con=XY(x,:);
        Con=[Con;WXY(i,:)];
        xc=sum(Con(:,1),1)/size(Con,1);
        yc=sum(Con(:,2),1)/size(Con,1);
        [t,R]=cart2pol(Con(:,1)-xc,Con(:,2)-yc);
        [P,n]=sort(t);
        Cn=Con(n,:);
        
        N=find((Cn(:,1)==WXY(i,1)).*(Cn(:,2)==WXY(i,2)));
        
        if N==1
            tet(i,1)=ygol(WXY(i,1),WXY(i,2),Cn(2,1),Cn(2,2), Cn(end,1), Cn(end,2));
        else
            if N==size(Cn,1)
                tet(i,1)=ygol(WXY(i,1),WXY(i,2),Cn(1,1),Cn(1,2), Cn(end-1,1), Cn(end-1,2));
            else
                tet(i,1)=ygol(WXY(i,1),WXY(i,2),Cn(N-1,1),Cn(N-1,2), Cn(N+1,1), Cn(N+1,2));
            end;
        end;
   %%%% на всякий случай     
%         if tet(i,1)<=pi/2
%             
%             T=ygol(WXY(i,1),WXY(i,2),Cn(1:end-1,1),Cn(1:end-1,2), Cn(2:end,1), Cn(2:end,2));
%             tet(i,1)=max(T);
%         end;
%         
%         if tet(i,1)<=pi/2
%             figure(1010); plot(Cn(:,1),Cn(:,2),'.');
%             hold on;
%             for k=1:size(Cn,1)
%                 text(Cn(k,1),Cn(k,2),num2str(k));
%             end;
%             plot(WXY(i,1),WXY(i,2),'O');
%         end;
    end;
    %%%%

    for j=1:size(x,1)
        l1(j)=L(x(j),y(j));
        b1(j)=B(x(j),y(j));
    end;

    Rk(i,1)=exp((sum(b1./l1.*log(l1),2)-tet(i,1))/sum(b1./l1,2));
    
    clear b1 l1 Con1
end;

for l=1:Nl
   Wcof(:,l)=K(Won,l).*H(Won,l)./log(Rk/r0).*Perf(:,l).*SDoly;
end;

Wf=2*pi*Wcof(:);


end


function TET=ygol(x1,y1,x2,y2,x3,y3)
xa=x2-x1;
ya=y2-y1;
xb=x3-x1;
yb=y3-y1;
TET=acos((xa.*xb+ya.*yb)./(((xa.^2+ya.^2).*(xb.^2+yb.^2)).^0.5));

end

