function Z=Gl_PRM
Nl=1;
Ns=5;                % Количество скважин
T=100;         % Период разработки.
dt=0;                % Шаг по времени % если 0- то динамический шаг
ndt=1;
g=9.81e-30;

%% ГДМ параметры
as=[2,2,2,2];                             % Фазовые в порах  (w-o) и (o-g)
aw=[0.1,1,1,1];                          % Крайняя точка фазовой в порах (w-o) и (o-g)
aw(3) = aw(2); 
ts=[1,1,1,1];                             % Фазовые в трещинах  (w-o) и (o-g)
tw=[1,1,1,10];                            % Крайняя точка фазовой в трещине (w-o) и (o-g)
kms=0*[1,1,1,1]*1e-3;                   % Коэф. для ур-я форгеймера поры-верт.тр. -гор.тр.-двойная ср.

SwrW=0.25;                               %Остаточная вода в системе (w-o), матрица
SorW=0.2625;                             %Остаточная нефть в системе (w-o), матрица

SorG=0.1;                               %Остаточная нефть в системе (o-g), матрица
SgrG=0;                             %Остаточный газ в системе (o-g), матрица

SwrWC=0;                                 %Остаточная вода в системе (w-o), трещина
SorWC=0;                                 %Остаточная нефть в системе (w-o),трещина

SorGC=0;                               %Остаточная нефть в системе (o-g), трещина
SgrGC=0;                             %Остаточный газ в системе (o-g), трещина

mu=[1,30,0.1,1];                       % Вязкости w-o-g-p
   %С      mu
mup=[0,  mu(1), 0.002,20];                  % вязкость полимера
Ro=[1.04,0.95,0.00066,3,2.7]*1000;      % Плотности вода-нефть-газ-песчаник-глина
zc=[1e-6,1e-4,1e-3,0,0,0]*0;   % Сжимаемость вода-нефть-газ-песчаник-глина-трещ            
Bo=[1,1,0.5,1,1];                       % Объёмные коэффициенты
rs = 0.01;

%% Параметры трещины
dh=0.05;                              % Раскрытость верт. дискретной трещины
hg=0.0005;                              % Раскрытость гор. дискретной трещины

Kc=1e+3;                                % Проницаемость дискретной трещины
Kg=1e+2;                                % Проницаемость горизонтальной трещины
Kd=1e+1;                                % Проницаемость двойной среды 

fC=[0,0,0];                             %Флаг на трещины верт.тр. - гор.тр. - двойная среда
ddol=0.1;				
drob=50;
Alp_C=100;                                % Коэф. обмена вертикальной трещины и пор

 %% Тепловые параметры 

lam=[0.65,0.13,0.06,0.5,0.8]*3600*24;    % Коэф. Теплопроводности
Cp=[4.2,1.9,2.22,0.75,0.75]*1000;        % Теплоёмкости


%% Запись в структуру
Z.Nl=Nl;
Z.as=as;
Z.aw=[aw,SwrW,SorW,SorG,SgrG];
Z.ts=ts;
Z.tw=[tw,SwrWC,SorWC,SorGC,SgrGC];
Z.mu=mu;
Z.Ns=Ns;
Z.Ta=T;
Z.dt=dt;
Z.ndt=ndt;
Z.Ro=Ro;
Z.lam=lam;
Z.Cp=Cp;
Z.dh=dh;
Z.Hg=hg;
Z.Kc=Kc;
Z.Kg=Kg;
Z.Kd=Kd;
Z.zc=zc;
Z.Bo=Bo;
Z.mup=[mup(1:2);mup(3:4)];
Z.kms=kms;
Z.g=g;
Z.fC=fC;
Z.drob=drob;
Z.ddol=ddol;
Z.Alp_C=Alp_C;
Z.rs=rs;

%[Z.Sc,Z.Fc]=sum2bolBO(Z.aw,Z.as,Z.mu);
%[Z.Sc2,Z.Fc2]=sum2bolBO(Z.tw,Z.ts,Z.mu);
[Z.Sc,Z.Fc]=sum2bol_Dima(Z.aw,Z.as,Z.mu);
[Z.Sc2,Z.Fc2]=sum2bol_Dima(Z.tw,Z.ts,Z.mu);
end
function [Sc,Fc]=sum2bol_Dima(AW,as,mu)
Swr=AW(5);
Sor=AW(6);
aw=AW(1:2);

Sw=Swr:0.001:(1-Sor);

SW=(Sw-Swr)./(1-Sor-Swr);
ko=(1-SW).^as(2);
kw=aw(1)*SW.^as(1);

gam=mu(1)/mu(2);

f=kw./(gam*ko+kw);
F=diff(f)./diff(Sw);
% figure(101);plot(Sw(2:end),F);
v=find(F==max(F));

Sc=Sw(v(1));
Fc=F(v(1));
end
function [Sc,Fc]=sum2bolBO(AW,as,mu)
SwrW=AW(5);
SorW=AW(6);
SorG=AW(7);
SgrG=AW(8);
aw=AW(1:4);

Sw=SwrW:0.001:(1-SorW);
Sg=SgrG:0.001:(1-SorG-SwrW);
SW=(Sw-SwrW)./(1-SorW-SwrW);
SG=(Sg-SgrG)./(1-SgrG-SorG-SwrW);
kg=aw(4)*SG.^as(4);
kog = aw(3)*(1-SG).^as(3);
kw=aw(1)*SW.^as(1);
kow = aw(2)*(1-SW).^as(2);

fw=kw./(kw+mu(1)/mu(2)*kow);
fg=kg./(kg+mu(3)/mu(2)*kog);
Fw=diff(fw)./diff(Sw);
Fg=diff(fg)./diff(Sg);
% figure(101);plot(Sw(2:end),F);
vw=find(Fw==max(Fw));
vg=find(Fg==max(Fg));

Sc=Sw(vw(1));
Fc=max(Fw(vw(1)),Fg(vg(1)));
end
