T=1:size(KIN,1);
t=0:365:T(end);

% KIN
% c
% sQo
% Vp
% sQz
% Qz
% sQl
% Ql

kin=interp1(T,KIN,t,'linear','extrap')';
C=interp1(T,c,t,'linear','extrap')';
sqo=interp1(T,sQo,t,'linear','extrap')';
vp=interp1(T,Vp,t,'linear','extrap')';
sqz=interp1(T,sQz,t,'linear','extrap')';
qz=interp1(T,Qz,t,'linear','extrap')';
sql=interp1(T,sQl,t,'linear','extrap')';
ql=interp1(T,Ql,t,'linear','extrap')';

T=T/365;