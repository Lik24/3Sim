function [A2CL,A2CW,A2CP,A2CG]=Obmen_T2M(A2C,Pm,Pt,Sw,Cw,K,PR,MCp,CCp)

as=PR.as;
aw=PR.aw;
ts=PR.ts;
tw=PR.tw;
mu=PR.mu;
mup=PR.mup;


[n,m]=size(A2C);
[r,c]=find(A2C);
A2C=A2C(:);
a2c=A2C(r+(c-1)*n).*K(r);

% Pm=P(1:n);
% Pt=P(n+1:end);

vP=Pm(r)>=Pt(c);
% vP
% [Pm(r),Pt(c)]
Swc=Sw(r);
Swl=Cw(c);

Cpc=MCp(r);
Cpl=CCp(c);

Swe=Swc.*vP+Swl.*(vP==0);
Cpe=Cpc.*vP+Cpl.*(vP==0);

wmu1=fusion(Cpe(vP==1),mup);
wmu0=fusion(Cpe(vP==0),mup);

Afw=Sat_cal(Swe(vP==1),1,1,as,aw)./wmu1; 
Afp=Sat_cal(Swe(vP==1),1,1,as,aw)./wmu1.*Cpe(vP==1); 
% Afw=Sat_cal(Swe(vP==1),1,1,as,aw).*((1-Cpe(vP==1))./mu(1)+Cpe(vP==1)./mu(4)); %water
Afo=Sat_cal(Swe(vP==1),2,1,as,aw)/mu(2); %oil
%Afp=Sat_cal(Swe(vP==1).*Cpe(vP==1),1,1,as,aw)/mu(4); %polim

Cfw=Sat_tube(Swe(vP==0),1,1,ts,tw)./wmu0; %water
Cfp=Sat_tube(Swe(vP==0),1,1,ts,tw)./wmu0.*Cpe(vP==0);
Cfo=Sat_tube(Swe(vP==0),2,1,ts,tw)/mu(2); %oil
%Cfp=Sat_tube(Swe(vP==0).*Cpe(vP==0),1,1,ts,tw)/mu(4); %polim

v1=find(vP==1);
v2=find(vP==0);
vone=ones(m,1);
% 
 %size([P])
% size([v1;v2])
% size(vone)
% size([Afw+Afo;Cfw+Cfo])

ACL=sparse([v1;v2],vone,[Afw+Afo;Cfw+Cfo],m,1);
ACW=sparse([v1;v2],vone,[Afw;Cfw],m,1);
ACP=sparse([v1;v2],vone,[Afp;Cfp],m,1);
%ACL
a2cl=a2c.*ACL;
a2cw=a2c.*ACW;
a2cp=a2c.*ACP;

A2CL=sparse(r,c,a2cl,n,m);
A2CW=sparse(r,c,a2cw,n,m);
A2CP=sparse(r,c,a2cp,n,m);
A2CG=1;