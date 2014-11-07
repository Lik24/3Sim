CD=load('history_deb_day');

T=1:size(CD.zak_lik2,2);
t=365/4:365/4:T(end);

Qo=sum(CD.dob_oil2,1);
Ql=sum(CD.dob_lik2,1);
Qz=sum(CD.zak_lik2,1);
C=1-Qo./Ql;

sQo=cumsum(Qo,1);
sQl=cumsum(Ql,1);
sQz=cumsum(Qz,1);

Nw_dob(:,1)=sum(CD.a0==1,1);
Nw_zak(:,1)=sum(CD.a0==-1,1);

Vp=sQz/sum(dV0);
V0=sum(dV0.*(1-Sw0));
KIN=sQo/V0;

PpW=CD.dob_Ppl2+CD.zak_Ppl2;

for t1=1:size(PpW,2) 
  Ppw(t1,1)=mean(PpW(PpW(:,t1)~=0,t1),1);
end
Pw=CD.dob_Pw2+CD.zak_Pw2;

W1=(CD.dob_lik2+CD.zak_lik2)./(PpW-Pw);
W1(isnan(W1)==1)=0;
for t1=1:size(Pw,2) 
    PW_d(t1,1)=mean(Pw(CD.a0(:,t1)==1,t1),1);
    PW_i(t1,1)=mean(Pw(CD.a0(:,t1)==-1,t1),1);
    W1_d(t1,1)=mean(W1(CD.a0(:,t1)==1,t1),1);
    W1_i(t1,1)=mean(W1(CD.a0(:,t1)==-1,t1),1);
end

% c(:,1)=interp1(T,C,t,'linear','extrap');    % - ������������� �������
% qo(:,1)=interp1(T,Qo,t,'linear','extrap');  % - ������� ������ �����,
% ql(:,1)=interp1(T,Ql,t,'linear','extrap');  % - ��������,
% qz(:,1)=-interp1(T,Qz,t,'linear','extrap'); % - �������
% sqo(:,1)=interp1(T,sQo,t,'linear','extrap')/1000; % - ����������� ������ �����, 
% sql(:,1)=interp1(T,sQl,t,'linear','extrap')/1000; % - ����������� ������ ��������,
% sqz(:,1)=interp1(T,sQz,t,'linear','extrap')/1000; % - ����������� �������
% 
% nw_dob(:,1)=interp1(T,Nw_dob,t,'linear','extrap');% - ����� ���. ���.,
% nw_zak(:,1)=interp1(T,Nw_zak,t,'linear','extrap');% - ����� ���. ���.
% 
% mqo=qo./nw_dob; % - ������� ����� �����
% mql=ql./nw_dob; % - ������� ����� ��������
% mqz=qz./nw_zak; % - ������� �����������
% soot=nw_dob./nw_zak; % - �����������
% skomp=sqz./sql; % - ����������� ���. � ���.
% tkomp=qz./ql;
% vp(:,1)=interp1(T,Vp,t,'linear','extrap');% - ����������� ������� �����
% kin(:,1)=interp1(T,KIN,t,'linear','extrap');% - ���
% ppw(:,1)=interp1(T,Ppw,t,'linear','extrap');% - ������� ��������� �������� � ����� �� �������, � ���� ������ (� �������, � ������� �������� ���������� ��������) � � ���� ������� - ���� ��� �� ������
% pwd(:,1)=interp1(T,PW_d,t,'linear','extrap');% - ������� �������� �������� �� ���������� ���.
% pwi(:,1)=interp1(T,PW_i,t,'linear','extrap');% - �������� �� ���.
% w1_d(:,1)=interp1(T,W1_d,t,'linear','extrap');% - ������� �������������� ���������� 
% w1_i(:,1)=interp1(T,W1_i,t,'linear','extrap');% - � �������������� �������

tt=[1:ceil(T(end)/365*4)]';
c(:,1)=osred(T,C',tt);    % - ������������� �������
qo(:,1)=osred(T,Qo',tt);  % - ������� ������ �����,
ql(:,1)=osred(T,Ql',tt);  % - ��������,
qz(:,1)=osred(T,Qz',tt); % - �������
sqo(:,1)=osred(T,sQo',tt)/1000; % - ����������� ������ �����, 
sql(:,1)=osred(T,sQl',tt)/1000; % - ����������� ������ ��������,
sqz(:,1)=osred(T,sQz',tt)/1000; % - ����������� �������

nw_dob(:,1)=osred(T,Nw_dob,tt);% - ����� ���. ���.,
nw_zak(:,1)=osred(T,Nw_zak,tt);% - ����� ���. ���.

mqo=qo./nw_dob; % - ������� ����� �����
mql=ql./nw_dob; % - ������� ����� ��������
mqz=qz./nw_zak; % - ������� �����������
soot=floor(nw_dob)./floor(nw_zak); % - �����������
skomp=sqz./sql; % - ����������� ���. � ���.
tkomp=qz./ql;
vp(:,1)=osred(T,Vp',tt);% - ����������� ������� �����
kin(:,1)=osred(T,KIN',tt);% - ���
ppw(:,1)=osred(T,Ppw,tt);% - ������� ��������� �������� � ����� �� �������, � ���� ������ (� �������, � ������� �������� ���������� ��������) � � ���� ������� - ���� ��� �� ������
pwd(:,1)=osred(T,PW_d,tt);% - ������� �������� �������� �� ���������� ���.
pwi(:,1)=osred(T,PW_i,tt);% - �������� �� ���.
w1_d(:,1)=osred(T,W1_d,tt);% - ������� �������������� ���������� 
w1_i(:,1)=osred(T,W1_i,tt);% - � �������������� �������

Text(1)={'�������'};
Text(2)={'������������� �������, �.��.'};
Text(3)={'������� ������ �����, �3/���'};
Text(4)={'������� ������ ��������, �3/���'};
Text(5)={'������� �������, �3/���'};

Text(6)={'�����. ������ �����, ���. �3'};
Text(7)={'�����. ������ ��������, ���. �3'};
Text(8)={'�����. �������, ���. �3'};

Text(9)={'����� ���. ���.'};
Text(10)={'����� ���. ���.'};
Text(11)={'�����������'};

Text(12)={'������� ����� �����, �3/���/���.'};
Text(13)={'������� ����� ��������, �3/���/���.'};
Text(14)={'������� �����������, �3/���/���.'};

Text(15)={'����������� ���., �.��.'};
Text(16)={'����������� ���., �.��.'};

Text(17)={'����������� ������� �����, ���.��.'};
Text(18)={'���, �.��.'};

Text(19)={'������� ��. ��������, ���.'};
Text(20)={'������� ���. �������� �� ���. ���., ���.'};
Text(21)={'������� ���. �������� �� ���. ���., ���.'};
Text(22)={'������� �������������� ���.'};
Text(23)={'������� �������������� ���'};
tt=1:size(c,1);
for i=1:size(kin,2)
    A(:,:,i)=[tt',c(:,i),qo(:,i),ql(:,i),qz(:,i),sqo(:,i),sql(:,i),sqz(:,i),floor(nw_dob(:,i)),floor(nw_zak(:,i)),soot(:,i),...
        mqo(:,i),mql(:,i),mqz(:,i),skomp(:,i),tkomp(:,i),vp(:,i),kin(:,i),ppw(:,i),...
        pwd(:,i),pwi(:,i),w1_d(:,i),w1_i(:,i)];
    xlswrite('Out_REZ_his.xls',Text,strcat('���_�',num2str(i)),'A1');
    xlswrite('Out_REZ_his.xls',A(:,:,i),strcat('���_�',num2str(i)),'A2');
end;