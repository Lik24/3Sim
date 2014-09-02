function [NT,PXY,gXY]=kvad_crack_fun(WXY,NL)
rad=40;         % ������ ������� ����� ����� �� ���������
drob=50;      % ��������� �����
dl=30;         %���������� ����� ������ �� �������
dl2=30;        %���������� ������� ����� ������ ������� ���� �������
ws=size(WXY,1);

for i=1:NL
    gg=[50,100;
        50,900];%�������������� �������
    %   gg=[50,950;450,550]; %��������� �������
    
    g_cr{1,1}=gg;
    [ p,GR ] = Mesh4( WXY,drob );
    [ crack2 ] = Fracture( g_cr,dl );
    [ pm2 ] = Purgatory( p,crack2,dl2,WXY,rad);
    cr=crack2{1,1};
    
   
    gXY=[cr;pm2];
%    nt=1:size(gXY,1);
    nt=1:size(cr,1);
    NT(i)={ws+nt};
    PXY(1,i)={[cr(1:end-1,1),cr(2:end,1)]};
    PXY(2,i)={[cr(1:end-1,2),cr(2:end,2)]};
end;


plot(pm2(:,1),pm2(:,2),'*');
hold on
plot(gg(:,1),gg(:,2));
plot(crack2{1,1}(:,1),crack2{1,1}(:,2),'o');

    % cr - ������ ����� �������
    % pm2 - ������ ����� ����� (��� �������)


