function [NLT,PXY,gXY,dl,tXY]=kvad_crack_fun(WXY,NL)
dh=2;

rad=6.5*dh;         % Радиус очистки вокруг скважины
drob=6.5*dh;        % Густота сетки
dl=6.5*dh;          % Шаг трещины
dl2=6.5*dh;         % Очитска вокруг трещины

ws=size(WXY,1);
NT=cell(NL,1);
PXY=cell(2,NL);
g_cr{1,1}=[];
   % Nt l    X  Y
g_cr{1,1}=[250,250;
           180,210]; %пїЅпїЅпїЅпїЅ

g_cr{2,1}=[210,210;
           140,170]; %пїЅпїЅпїЅпїЅ
       
%g_cr{2,2}=[110,110;
%           140,170]; %пїЅпїЅпїЅпїЅ

ImA=zeros(size(g_cr,1),NL);

if isempty(g_cr{1,1})==0; 
    [p,GR] = Mesh4(WXY,drob);

    for i=1:size(g_cr,1)
        for j=1:size(g_cr,2)
            A=g_cr{i,j};
            ImA(i,j)=(isempty(A)==0);
        end;
    end;

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

    gXY=[cr;pm2];
    gXY=NODuble2(gXY);
    
    sr=[];
    for j=1:ws
        rr=find((WXY(j,1)==gXY(:,1)).*(WXY(j,2)==gXY(:,2))==1);
        sr=[sr;rr];
    end;
    
    if isempty(sr)==0
        gXY(sr,:)=[];
    end;
    tXY=[WXY;gXY];
        

for l=1:NL
rl=find(ImA(:,l));  
  for i=rl 
     cr=crk{i,l};
     nt=[];
      for j=1:size(cr,1)
       a1=(cr(j,1)==tXY(:,1)).*(cr(j,2)==tXY(:,2));
       nt(j)=find(a1==1);
      end;
     Nt(i)={nt};
  end;
NLT(l)={Nt};
Nt={[]};
end;

    
for l=1:NL
rl=find(ImA(:,l));  
  for i=rl 
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