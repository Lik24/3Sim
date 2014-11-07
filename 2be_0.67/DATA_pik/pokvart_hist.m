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

% c(:,1)=interp1(T,C,t,'linear','extrap');    % - обводненность текущая
% qo(:,1)=interp1(T,Qo,t,'linear','extrap');  % - Текущая добыча нефти,
% ql(:,1)=interp1(T,Ql,t,'linear','extrap');  % - Жидкости,
% qz(:,1)=-interp1(T,Qz,t,'linear','extrap'); % - Закачка
% sqo(:,1)=interp1(T,sQo,t,'linear','extrap')/1000; % - Накопленные добыча нефти, 
% sql(:,1)=interp1(T,sQl,t,'linear','extrap')/1000; % - Накопленные добыча жидкости,
% sqz(:,1)=interp1(T,sQz,t,'linear','extrap')/1000; % - Накопленные закачка
% 
% nw_dob(:,1)=interp1(T,Nw_dob,t,'linear','extrap');% - Число доб. скв.,
% nw_zak(:,1)=interp1(T,Nw_zak,t,'linear','extrap');% - Число наг. скв.
% 
% mqo=qo./nw_dob; % - Средний дебит нефти
% mql=ql./nw_dob; % - Средний дебит жидкости
% mqz=qz./nw_zak; % - Средняя приёмистость
% soot=nw_dob./nw_zak; % - соотношение
% skomp=sqz./sql; % - Компенсация тек. и нак.
% tkomp=qz./ql;
% vp(:,1)=interp1(T,Vp,t,'linear','extrap');% - Прокаченный поровый объем
% kin(:,1)=interp1(T,KIN,t,'linear','extrap');% - КИН
% ppw(:,1)=interp1(T,Ppw,t,'linear','extrap');% - Среднее пластовое давление в целом по участку, в зоне отбора (в ячейках, в которые попадают добывающие скважины) и в зоне закачки - если это не сложно
% pwd(:,1)=interp1(T,PW_d,t,'linear','extrap');% - Среднее забойное давление по добывающим скв.
% pwi(:,1)=interp1(T,PW_i,t,'linear','extrap');% - отдельно по наг.
% w1_d(:,1)=interp1(T,W1_d,t,'linear','extrap');% - средняя продуктивность добывающих 
% w1_i(:,1)=interp1(T,W1_i,t,'linear','extrap');% - и нагнетательных скважин

tt=[1:ceil(T(end)/365*4)]';
c(:,1)=osred(T,C',tt);    % - обводненность текущая
qo(:,1)=osred(T,Qo',tt);  % - Текущая добыча нефти,
ql(:,1)=osred(T,Ql',tt);  % - Жидкости,
qz(:,1)=osred(T,Qz',tt); % - Закачка
sqo(:,1)=osred(T,sQo',tt)/1000; % - Накопленные добыча нефти, 
sql(:,1)=osred(T,sQl',tt)/1000; % - Накопленные добыча жидкости,
sqz(:,1)=osred(T,sQz',tt)/1000; % - Накопленные закачка

nw_dob(:,1)=osred(T,Nw_dob,tt);% - Число доб. скв.,
nw_zak(:,1)=osred(T,Nw_zak,tt);% - Число наг. скв.

mqo=qo./nw_dob; % - Средний дебит нефти
mql=ql./nw_dob; % - Средний дебит жидкости
mqz=qz./nw_zak; % - Средняя приёмистость
soot=floor(nw_dob)./floor(nw_zak); % - соотношение
skomp=sqz./sql; % - Компенсация тек. и нак.
tkomp=qz./ql;
vp(:,1)=osred(T,Vp',tt);% - Прокаченный поровый объем
kin(:,1)=osred(T,KIN',tt);% - КИН
ppw(:,1)=osred(T,Ppw,tt);% - Среднее пластовое давление в целом по участку, в зоне отбора (в ячейках, в которые попадают добывающие скважины) и в зоне закачки - если это не сложно
pwd(:,1)=osred(T,PW_d,tt);% - Среднее забойное давление по добывающим скв.
pwi(:,1)=osred(T,PW_i,tt);% - отдельно по наг.
w1_d(:,1)=osred(T,W1_d,tt);% - средняя продуктивность добывающих 
w1_i(:,1)=osred(T,W1_i,tt);% - и нагнетательных скважин

Text(1)={'Квартал'};
Text(2)={'Обводненность текущая, д.ед.'};
Text(3)={'Текущая добыча нефти, м3/сут'};
Text(4)={'Текущая добыча жидкости, м3/сут'};
Text(5)={'Текущая закачка, м3/сут'};

Text(6)={'Накоп. добыча нефти, тыс. м3'};
Text(7)={'Накоп. добыча жидкости, тыс. м3'};
Text(8)={'Накоп. закачка, тыс. м3'};

Text(9)={'Число доб. скв.'};
Text(10)={'Число наг. скв.'};
Text(11)={'Соотношение'};

Text(12)={'Средний дебит нефти, м3/сут/скв.'};
Text(13)={'Средний дебит жидкости, м3/сут/скв.'};
Text(14)={'Средняя приёмистость, м3/сут/скв.'};

Text(15)={'Компенсация нак., д.ед.'};
Text(16)={'Компенсация тек., д.ед.'};

Text(17)={'Прокаченный поровый объем, отн.ед.'};
Text(18)={'КИН, д.ед.'};

Text(19)={'Среднее пл. давление, атм.'};
Text(20)={'Среднее заб. давление по доб. скв., атм.'};
Text(21)={'Среднее заб. давление по наг. скв., атм.'};
Text(22)={'Средняя продуктивность доб.'};
Text(23)={'Средняя продуктивность наг'};
tt=1:size(c,1);
for i=1:size(kin,2)
    A(:,:,i)=[tt',c(:,i),qo(:,i),ql(:,i),qz(:,i),sqo(:,i),sql(:,i),sqz(:,i),floor(nw_dob(:,i)),floor(nw_zak(:,i)),soot(:,i),...
        mqo(:,i),mql(:,i),mqz(:,i),skomp(:,i),tkomp(:,i),vp(:,i),kin(:,i),ppw(:,i),...
        pwd(:,i),pwi(:,i),w1_d(:,i),w1_i(:,i)];
    xlswrite('Out_REZ_his.xls',Text,strcat('Вар_№',num2str(i)),'A1');
    xlswrite('Out_REZ_his.xls',A(:,:,i),strcat('Вар_№',num2str(i)),'A2');
end;