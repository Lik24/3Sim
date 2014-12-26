function [A2CW,A2CP,A2CG]=Obmen_T2M_1(a2c,P,Sw,Cw,ts,tw,as,aw,mu,MCp,CCp,r,c)
n=size(Sw,1);
m=size(Cw,1);

Pm=P(1:n);
Pt=P(n+1:end);

vP=Pm(r)>=Pt(c);
Swc=Sw(r);
Swl=Cw(c);

Cpc=MCp(r);
Cpl=CCp(c);

Swe=Swc.*vP+Swl.*(vP==0);
Cpe=Cpc.*vP+Cpl.*(vP==0);

Af=Sat_cal(Swe(vP==1),1,1,as,aw); %water
Cf=Sat_tube(Swe(vP==0),1,1,ts,tw); %water

Afw=Af.*((1-Cpe(vP==1))./mu(1)+Cpe(vP==1)/mu(4));
Afp=Af.*Cpe(vP==1)/mu(4);

Cfw=Cf.*((1-Cpe(vP==0))./mu(1)+Cpe(vP==0)/mu(4));
Cfp=Cf.*Cpe(vP==0)./mu(4);

v1=find(vP==1);
v2=find(vP==0);
vone=ones(m,1);

ACW=sparse([v1;v2],vone,[Afw;Cfw],m,1);
ACP=sparse([v1;v2],vone,[Afp;Cfp],m,1);

% ACW=op_accuamarray([v1;v2],m);
% ACW=ACW.*[Afw+Afp;Cfw+Cfp];
% ACP=op_accuamarray([v1;v2],[Afp;Cfp]);

A2CW=a2c.*ACW;
A2CP=a2c.*ACP;

A2CG=1;