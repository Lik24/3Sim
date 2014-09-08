function Z=Gl_PRM
Nl=1;
Ns=4;                             % Количество скважин
T=365*3;              % Период разработки.
dt=0.5;                % Шаг по времени % если 0- то динамический шаг
Ta=T;
ndt=1;

%% ГДМ параметры
as=[2,2,1];                             % Фазовые в порах  w-o-g
aw=[0.1,1,10];                          % Крайняя точка фазовой в порах w-o-g
 
ts=[2,2,1];                             % Фазовые в трещинах  w-o-g
tw=[0.1,1,10];                            % Крайняя точка фазовой в трещине w-o-g
bet=0*1e-2;                             % Коэф. для ур-я форгеймера

mu=[1,300,0.1,1];                       % Вязкости w-o-g-p
Ro=[1.04,0.95,0.00066,3,2.7]*1000;      % Плотности вода-нефть-газ-песчаник-глина
zc=[1e-6,1e-5,1e-3,1e-6,1e-6];          % Сжимаемость             
Bo=[1,0.95,0.5,1,1];                    % Объёмные коэффициенты

%% Параметры трещины
dh=0.0005;
Kc=1e+3;                                % Проницаемость трещины

%% Тепловые параметры 
lam=[0.65,0.13,0.06,0.5,0.8]*3600*24;    % Коэф. Теплопроводности
Cp=[4.2,1.9,2.22,0.75,0.75]*1000;        % Теплоёмкости
%Cp=[4.2,4.2,2.22,0.75,0.75]*1000;       % Теплоёмкости


%% Запись в структуру
Z.Nl=Nl;
Z.as=as;
Z.aw=aw;
Z.ts=ts;
Z.tw=tw;
Z.mu=mu;
Z.Ns=Ns;
Z.Ta=Ta;
Z.dt=dt;
Z.ndt=ndt;
Z.Ro=Ro;
Z.bet=bet;
Z.lam=lam;
Z.Cp=Cp;
Z.dh=dh;
Z.Kc=Kc;
Z.zc=zc;
Z.Bo=Bo;

%% Символьный
syms Sw
ko=(1-Sw).^as(2);
kw=aw(1)*Sw.^as(1);

gam=mu(1)/mu(2);

f=kw/(gam*ko+kw);
F=diff(f,Sw);

dS=0:0.001:1;
fn=(subs(f,dS)-subs(f,0))./(dS-0);

fm=max(eval(fn));%
fn=eval(fn);
[r,c]=find(fm==fn);

Sc=dS(c);
Fc=eval(subs(F,Sc));

kot=(1-Sw).^ts(2);
kwt=tw(1)*Sw.^ts(1);

f=kwt/(gam*kot+kwt);
F=diff(f,Sw);

ds=0.001;
dS=0:ds:1;
fn=(subs(f,dS)-subs(f,0))./(dS-0);
fm=max(eval(fn));%
fn=eval(fn);

[r,c]=find(fm==fn);
Sc2=dS(c(1));

if Sc2==ds
    Sc2=mu(1)/mu(2)*tw(2)/tw(1);
    Fc2=eval(subs(F,Sc2));
else
    Fc2=eval(subs(F,Sc2));
end;

Z.Fc=Fc;
Z.Fc2=Fc2;
Z.Sc2=Sc2;
