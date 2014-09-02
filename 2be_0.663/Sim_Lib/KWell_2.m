function Wf=KWell_2(K,H,S,L,B,Won,r,c,Dolya,r0)
Perf=1;

tet=Dolya.*2.*pi;
 %tet=[1;0.25;0.25;0.25;0.25].*2.*pi; %для пятиточки
for i=1:size(Won,1);
    x=r(c==Won(i));
    y=c(c==Won(i));
        for j=1:size(x,1)
              l1(j)=L(x(j),y(j));
              b1(j)=B(x(j),y(j));
        end;
%     sum(b1./l1.*log(l1),2)
%     tet(i,1)
%     sum(b1./l1,2)
% exp((sum(b1./l1.*log(l1),2)-tet(i,1))/sum(b1./l1,2))

Rk(i,1)=exp((sum(b1./l1.*log(l1),2)-tet(i,1))/sum(b1./l1,2));

clear b1 l1
end;


Wcof=K(Won).*H(Won)./log(Rk/r0)*Perf;
Wf=2*pi*Wcof;
