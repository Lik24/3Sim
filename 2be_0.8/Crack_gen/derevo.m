function [nt,PXY]=derevo(nt,gXY,ii)

PXY=cell(2,size(nt,1));

for l=1:size(nt,1)
    [A]=MR_Prop(gXY,1);
    na=size(gXY,1);
    Nt=nt{l}-(l-1)*na;
    new_nt=Nt;
    SA=1:size(A,1);
  
    A(:,Nt)=0;
%     x1=zeros(1,2);
%     y1=zeros(1,2);
%     x1(1,:)=[];
%     y1(1,:)=[];

% size(gXY(Nt(1:end-1),1))
% size(gXY(Nt(2:end),1))

    x1=[gXY(Nt(1:end-1),1),gXY(Nt(2:end),1)];
    y1=[gXY(Nt(1:end-1),2),gXY(Nt(2:end),2)];
     figure(ii+l),plot(gXY(Nt,1),gXY(Nt,2),0,0,250,250);
     hold on
     jj=0;
     for i=2:size(Nt,2)-1
         %A=A1;
         
         tc=Nt(i);
         fl=0;
         j=0;
         
         while fl==0
             j=j+1;
             jj=jj+1;
             tc1=tc(1);
             [tc,nj,A]=nos_nof(tc(1),SA,A);
             
             new(j:j+size(nj,2)-1)=nj;
             
             j=j+size(nj)-1;
             fl=j>size(Nt,2)-i+10;
             
             tc1=repmat(tc1,size(tc,2),1);
             figure(ii+l),plot([gXY(tc1,1)';gXY(tc,1)'],[gXY(tc1,2)';gXY(tc,2)'],'b');
            
             x1=[x1;[gXY(tc1,1),gXY(tc,1)]];
             y1=[y1;[gXY(tc1,2),gXY(tc,2)]];
             
             hold on
         end
         
         
         new(new==-1)=[];
         new_nt=[new_nt,new];
       
     end;
new_nt=NODuble(new_nt);
nt(l)={new_nt+(l-1)*na};

PXY(1,l)={x1};
PXY(2,l)={y1};

end;

end

function [tc,new,A]=nos_nof(tc,SA,A)
cn=SA(A(tc,:)==1);
if isempty(cn)==0
    nn=randi(1:2,1);
    nb=randi(size(cn,2),1,nn);
    tc=cn(nb);
    new=tc;
    A(:,tc)=0;
else
    new=-1;
    %A(:,tc)=0;
end;
end

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