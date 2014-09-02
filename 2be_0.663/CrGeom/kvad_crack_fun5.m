function [NT,PXY,gXY,dl,tXY]=kvad_crack_fun5(WXY,NL)
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
% g_cr{1,1}=[mean(WXY(:,1))-50,mean(WXY(:,2))+50;
%            mean(WXY(:,1))-70,mean(WXY(:,2))+70];%пїЅпїЅпїЅ\

g_cr{1,1}=[50,200;
           200,50]; %пїЅпїЅпїЅпїЅ

%g_cr{2,1}=[750,750;
%           250,250]; %пїЅпїЅпїЅпїЅ

 %g_cr{2,2}=[250,750;
 %           750,250]; %пїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅпїЅ пїЅпїЅпїЅпїЅпїЅпїЅпїЅ

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
        [pm2] = Purgatory( p,crk2,dl2,WXY,rad);

    cr=[];    
    for i=r'
        ds=crk2{i,1};
        cr=[cr;ds];
    end;
    gXY=[cr;pm2];
    gXY=NODuble2(gXY);
    
    for i=1:max(c)
        r1=r(c==i);
        cr=[];
        x1=[];
        x2=[];
        y1=[];
        y2=[];

        for j=r1'
            crn=crk2{j,1};
            cr=[cr;crn];

            x1=[x1;crn(1:end-1,1)];
            x2=[x2;crn(2:end,1)];
            y1=[y1;crn(1:end-1,2)];
            y2=[y2;crn(2:end,2)];
        end;

        nt=1:size(cr,1);
        sr=[];
        nj=[];
        for j=1:ws
            rr=find((WXY(j,1)==gXY(:,1)).*(WXY(j,2)==gXY(:,2))==1);
            sr=[sr;rr];
            nj=[nj;j*ones(size(rr,1),1)];
        end;

        if isempty(sr)==0
            gXY(sr,:)=[];
%             for j=1:size(sr,1)
%                 rr=find(sr(j)==nt);
%                 if isempty(rr)==0
%                     nt(rr+1:end)=nt(rr+1:end)-1;
%                     nt(rr)=nj(j)-ws;
%                     %nt(sr)=nj-ws;
%                 end;
%             end;
        end;
        
        tXY=[WXY;gXY];
        nc=size(tXY,1);
        
        for j=1:size(cr,1)
            a1=(cr(j,1)==tXY(:,1)).*(cr(j,2)==tXY(:,2));
            fg=find(a1==1);
            nt(j)=fg;
        end;
        nt=NODuble(nt);
 
        NT(i)={nc*(i-1)+nt};
        PXY(1,i)={[x1,x2]};
        PXY(2,i)={[y1,y2]};

    end;


    for i=1:size(gcr,1)
        gg=gcr{i,1};
    %     plot(gg(:,1),gg(:,2));
    %    plot(crk2{i,1}(:,1),crk2{i,1}(:,2),'o');
    end;
else
       [pm2] = Purgatory( p,crk2,dl2,WXY,rad);
end;
   
%     for i=1:2
%         plot(tXY(NT{i}-nc*(i-1),1),tXY(NT{i}-nc*(i-1),2))
%         hold on
%     end;
%     nc
%  rtyuj 

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