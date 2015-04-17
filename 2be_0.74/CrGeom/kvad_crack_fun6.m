function [NT,PXY,gXY]=kvad_crack_fun5(WXY,NL)
rad=6.5*2*2;         % ������ ������� ����� ����� �� ���������
drob=6.5*2*2;      % ��������� �����
dl=6.5*2*2;         %���������� ����� ������ �� �������
dl2=6.5*2;        %���������� ������� ����� ������ ������� ���� �������
ws=size(WXY,1);
NT=cell(NL,1);
g_cr{1,1}=[];
  % Nt l    X  Y
g_cr{1,1}=[250,750;
           750,250];%���\

%g_cr{2,1}=[250,250;
%          500,500]; %����

% g_cr{3,2}=[450,950;
%            450,50]; %��������� �������

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
    nc=size(gXY,1);

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
            for j=1:size(sr,1)
                rr=find(sr(j)==nt);
                if isempty(rr)==0
                    nt(rr+1:end)=nt(rr+1:end)-1;
                    nt(rr)=nj(j)-ws;
                    %nt(sr)=nj-ws;
                end;
            end;
        end;

        nt=NODuble(nt);

% nc*(i-1)+ws+nt
% nc*(i-1)
% ws
% nt
%gXY(nc*(i-1)+ws+nt,:)
%jkhkjh
        NT(i)={nc*(i-1)+ws+nt};
        PXY(1,i)={[x1,x2]};
        PXY(2,i)={[y1,y2]};

    end;

    % plot(pm2(:,1),pm2(:,2),'*');
    % hold on

    for i=1:size(gcr,1)
        gg=gcr{i,1};
    %     plot(gg(:,1),gg(:,2));
    %    plot(crk2{i,1}(:,1),crk2{i,1}(:,2),'o');
    end;
    else
        [pm2] = Purgatory( p,crk2,dl2,WXY,rad);
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
