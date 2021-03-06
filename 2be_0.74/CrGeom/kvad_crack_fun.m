function [Nt,gXY,dl,tXY,XY_GY,Won,WData]=kvad_crack_fun(XY_GY,NL,WData,dh,Won)

WXY=WData.WXY;
Won_s=Won;

rad=dh;         % ������ ������� ������ ��������
drob=dh;        % ������� �����
dl=dh;          % ��� �������
dl2=dh;         % ������� ������ �������
kx=1;
l=100;

g_cr{1,1}=[];
   % Nt l    X  Y  Z
   % ��������� ������������ �������
%    si=20;
%      kx=1+(rand(si,1)-0.5);
%      ax1=randi(max(XY_GY(:,1))/20,si,1)*20;
%      ay1=randi(max(XY_GY(:,2))/20,si,1)*20;
%      b1=ay1-kx.*ax1;
%      a=kx.^2+1;
%      b=-2*ax1-2*ay1.*kx+2.*kx.*b1;
%      c=-(l^2-ax1.^2-ay1.^2+2*ay1.*b1-b1.^2);
%      d=b.^2-4*a.*c;
%      ax2=(-b+d.^0.5)/2./a;
%      ay2=kx.*ax2+b1;
%               
%      [axy1,ia,~]=unique([ax1,ay1],'rows');
%      ax2=ax2(ia); ay2=ay2(ia);
%      r=[];
%      for i=1:size(ax2,1)
%        [k1,d1]=dsearchn(axy1(i,:),[ax2,ay2]);
%        r=[r;find(d1<dh*2)];
%      end
%      axy1(r,:)=[]; ax2(r)=[];   ay2(r)=[];
%      in=inpolygon(ax2,ay2,XY_GY(:,1),XY_GY(:,2));
%      
%   for i=1:size(ax2,1)  
%      g_cr{i}=[axy1(i,:),1;
%               ax2(i),ay2(i),1];  
%           plot([axy1(i,1),ax2(i)],[axy1(i,2),ay2(i)])
%           hold on
%           
%    end
 %  g_cr(in==0)=[];
 
 g_cr{1}=[0,0,1
          200,200,1]; %����
% 
%  g_cr{2}=[450,425,1;
%           280,305,1]; %����
% % %        
%  g_cr{3}=[395,350,1;
%           200,225,1]; %����
% %        
%  g_cr{4}=[120,125,1;
%           250,250,1]; %����
% 
% g_cr{5,1}=[100,510;
%            900,510]; %����
% %        
% g_cr{6,1}=[100,550;
%            900,550]; %����
%        
% g_cr{7,1}=[100,590;
%            900,590]; %����
% %        
% g_cr{8,1}=[100,630;
%            900,630]; %����
%        
% g_cr{9,1}=[100,670;
%            900,670]; %����
%        
% g_cr{10,1}=[100,710;
%             900,710]; %����
%         
% g_cr{11,1}=[100,755;
%            900,755]; %����
%        
% g_cr{12,1}=[100,790;
%            900,790]; %����
%        
% g_cr{13,1}=[100,830;
%             900,830]; %����

if isempty(g_cr{1,1})==0;
    [p,~] = Mesh3(XY_GY,drob);
    [WXY,Won,gcr,de]=crack2horwell(g_cr,WXY,Won,drob);
    
    [crk2] = Fracture(gcr,dl);
    [pm2] = Purgatory(p,crk2,dl2,WXY,rad,1,1);
    cr=cell2mat(crk2);
    
    [IN,ON]=inpolygon(pm2(:,1),pm2(:,2),XY_GY(:,1),XY_GY(:,2));
    gXY=[cr;pm2;pm2(ON==1,:)];
    gXY=unique(gXY,'rows');
    XY_GY=p(ON==1,:);
    
    [~,~,cb]=intersect(WXY,gXY,'rows');
    gXY(cb,:)=[];
    tXY=[WXY;gXY];
    n=size(tXY,1);

    for i=1:size(crk2,1)
        cr=crk2{i};
        % nt=[];
        [~,~,nt]=intersect(cr,tXY,'rows','stable');
        nt1=nt(1:end-1)';
        nt2=nt(2:end)';
        Z=gcr{i}(:,3);
        z(1)=min(Z).*(min(Z)<=NL)+NL*(min(Z)>NL);
        z(2)=max(Z).*(max(Z)<=NL)+NL*(max(Z)>NL);
        Z_nt1=[];
        Z_nt2=[];
        nl=size(nt1,2);
        k=0;
        for j=z(1):z(2)
            k=k+1;
            Z_nt1=[Z_nt1,nt1+(j-1)*n];
            Z_nt2=[Z_nt2,nt2+(j-1)*n];
        end
        if z(1)~=z(2)
            Z_nt1=[Z_nt1,Z_nt1(1:nl*(k-1))];
            Z_nt2=[Z_nt2,Z_nt1(nl+1:nl*k)];
        end
        Nt(i)={[Z_nt1;Z_nt2]};
    end;

else
    
    [p,GR] = Mesh3(XY_GY,drob);
    [WXY,Won,gcr,de]=crack2horwell(g_cr,WXY,Won,drob);
    [p] = Purgatory(p,0,dl2,WXY,rad,0,1);
    
    gXY=p;
    tXY=[WXY;gXY];
    Nt={zeros(2,0)};
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