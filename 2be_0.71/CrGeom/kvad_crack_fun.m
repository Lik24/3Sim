function [NLT,PXY,gXY,dl,tXY,XY_GY,Won,WData]=kvad_crack_fun(XY_GY,NL,WData,dh,Won)

WXY=WData.WXY;
Won_s=Won;

rad=dh;         % Радиус очистки вокруг скважины
drob=dh;        % Густота сетки
dl=dh;          % Шаг трещины
dl2=dh;         % Очитска вокруг трещины


NT=cell(NL,1);
PXY=cell(NL,1);

g_cr{1,1}=[];
   % Nt l    X  Y
g_cr{1,1}=[100,350;
           400,350]; %пїЅпїЅпїЅпїЅ

% g_cr{2,1}=[100,390;
%            900,390]; %пїЅпїЅпїЅпїЅ
% %        
% g_cr{3,1}=[100,430;
%            900,430]; %пїЅпїЅпїЅпїЅ
%        
% g_cr{4,1}=[100,470;
%            900,470]; %пїЅпїЅпїЅпїЅ
% 
% g_cr{5,1}=[100,510;
%            900,510]; %пїЅпїЅпїЅпїЅ
% %        
% g_cr{6,1}=[100,550;
%            900,550]; %пїЅпїЅпїЅпїЅ
%        
% g_cr{7,1}=[100,590;
%            900,590]; %пїЅпїЅпїЅпїЅ
% %        
% g_cr{8,1}=[100,630;
%            900,630]; %пїЅпїЅпїЅпїЅ
%        
% g_cr{9,1}=[100,670;
%            900,670]; %пїЅпїЅпїЅпїЅ
%        
% g_cr{10,1}=[100,710;
%             900,710]; %пїЅпїЅпїЅпїЅ
%         
% g_cr{11,1}=[100,755;
%            900,755]; %пїЅпїЅпїЅпїЅ
%        
% g_cr{12,1}=[100,790;
%            900,790]; %пїЅпїЅпїЅпїЅ
%        
% g_cr{13,1}=[100,830;
%             900,830]; %пїЅпїЅпїЅпїЅ

ImA=zeros(size(g_cr,1),NL);

if isempty(g_cr{1,1})==0; 
        
    [p,GR] = Mesh3(XY_GY,drob);
%plot(p(:,1),p(:,2),'*')
    for i=1:size(g_cr,1)
        for j=1:size(g_cr,2)
            A=g_cr{i,j};
            ImA(i,j)=(isempty(A)==0);
        end;
    end;

    [WXY,Won,g_cr,de]=crack2horwell(g_cr,WXY,Won,drob);
    ws=size(WXY,1);
    [r,c]=find(ImA);
    for i=1:size(r,1)
        gcr(i,1)=g_cr(r(i),c(i));
    end;

     
    [crk2] = Fracture(gcr,dl);
    [pm2] = Purgatory(p,crk2,dl2,WXY,rad);
    cr=cell2mat(crk2);

    for i=1:size(r,1)
       crk(r(i),c(i))=crk2(i);
    end;

    [IN,ON]=inpolygon(pm2(:,1),pm2(:,2),XY_GY(:,1),XY_GY(:,2));
    gXY=[cr;pm2;pm2(ON==1,:)];
    gXY=unique(gXY,'rows');
    XY_GY=p(ON==1,:);
    
    [~,~,cb]=intersect(WXY,gXY,'rows');
    gXY(cb,:)=[];
    tXY=[WXY;gXY];
%    tXY=unique(tXY,'rows');
%    [WXY1,ca,cb]=intersect(WXY,tXY,'rows','stable');

for l=1:NL
rl=find(ImA(:,l)); 
  for i=rl' 
     cr=crk{i,l};
    % nt=[];
     [~,~,nt]=intersect(cr,tXY,'rows','stable');
     nt1=nt(1:end-1)';
     nt2=nt(2:end)';
     Nt(i)={[nt1;nt2]};
  end;
NLT(l)={Nt};
Nt={[]};
end;

   
for l=1:NL
rl=find(ImA(:,l));  
  for i=rl' 
   crn=crk{i,l};
   x1=crn(1:end-1,1);
   x2=crn(2:end,1);
   y1=crn(1:end-1,2);
   y2=crn(2:end,2);
            
   pxy(1,i)={[x1,x2]};
   pxy(2,i)={[y1,y2]};
  end;
PXY(l)={pxy};
pxy={[]};
end;

else
    
   [p,GR] = Mesh4(WXY,drob); 
   NLT=cell(Nl,1);
   PXY=[];
   gXY=p;
   tXY=[WXY;gXY];
end;
   
doly=WData.Doly;
Sdoly=WData.SDoly;

won_s=Won_s(Won_s(:,2)==1,:);
doly_s=doly(Won_s(:,2)==1,:);
Sdoly_s=Sdoly(Won_s(:,2)==1);
do1y=doly;
Sdo1y=Sdoly;
do1y(de,:)=[];
Sdo1y(de,:)=[];

wn=unique(won_s(:,1));
 for i=1:size(wn,1)
     df=won_s(won_s(:,1)==wn(i),3);
     nf=Won(Won(:,1)==wn(i),3);
     do1y(Won(:,1)==wn(i),:)=interp1(df,doly_s(won_s(:,1)==wn(i),:),nf,'linear','extrap');
     Sdo1y(Won(:,1)==wn(i))=interp1(df,Sdoly_s(won_s(:,1)==wn(i)),nf,'linear','extrap');
 end

WData.WXY=WXY;
WData.Doly=do1y;
WData.SDoly=Sdo1y;
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

function XY=NODuble2(XY)
i=0;
while i<size(XY,1)
    i=i+1;
    a1=(XY(i,1)==XY(:,1)).*(XY(i,2)==XY(:,2));
    if sum(a1)>1
     r=find(a1);
     XY(r(2:end),:)=[];
    end;
end;
end