function [nt,PXY]=gls_crk(nt,gXY,kol,dlin,wfl,ii)

PXY=cell(2,size(nt,1));

for l=1:size(nt,1)
    [A]=MR_Prop(gXY,1);
    na=size(gXY,1);
    Nt=nt{l}-(l-1)*na;
    new_nt=[];
    SA=1:size(A,1);
    A=A.*(A>0);
    %A(:,Nt)=0;
x3=[];
y3=[];
    x1=gXY(:,1);
    y1=gXY(:,2);
     figure(ii+l),plot(min(gXY(:,1)),min(gXY(:,2)),max(gXY(:,1)),max(gXY(:,2)));
     hold on
 
     for i=1:kol
        k=0; 
        fl1=1;
        while fl1==1 && k<20 
         k=k+1;   
         ni=randi(na,1);
          if sum(ni==new_nt)==0
            fl1=0;  
          end;
        end;
        
         new=ni;
         for j=1:dlin
             r=find(A(:,ni));
             if isempty(r)==0
                 kv=randi(3,1);
                 % mk=min([max(kv),numel(r)]);
                 xr=rand(1);
                 mk=2*(xr>0.6)+(xr<=0.6);
                 nj=randi(numel(r),mk,1);
                 A(r(nj),:)=0;
                 
                 x2=x1(ni)*ones(1,mk);
                 y2=y1(ni)*ones(1,mk);
                 
                 
                 plot([x2;x1(r(nj))'],[y2;y1(r(nj))']);
                 ni=r(nj(1));

                 x3=[x3;[x2',x1(r(nj))]];
                 y3=[y3;[y2',y1(r(nj))]];
               
                 new=[new,r(nj)'];
             end;
         end;
          new_nt=[new_nt,new];  
     end;
     
 %new_nt=[new_nt,new];    
 new_nt=NODuble(new_nt);
 nt(l)={new_nt+(l-1)*na};    
     
 PXY(1,l)={x3};
 PXY(2,l)={y3};   
     
 
 end;
 hold off
 end
% 
% function [tc,new,A]=nos_nof(tc,SA,A)
% cn=SA(A(tc,:)==1);
% if isempty(cn)==0
%     nn=randi(1:2,1);
%     nb=randi(size(cn,2),1,nn);
%     tc=cn(nb);
%     new=tc;
%     A(:,tc)=0;
% else
%     new=-1;
%     %A(:,tc)=0;
% end;
% end
% 
function XY=NODuble(XY)
i=0;
while i<size(XY,2)
    i=i+1;
    a1=XY(i)==XY;
    if sum(a1)>1
     r=find(a1);
     XY(r(2:end))=[];
    end;
end;
end