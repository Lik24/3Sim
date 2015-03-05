function [A2CL,A2CW,A2CO,A2CG,A2CWCgs,A2COCgs,A2CP]=Obmen_T2M(A2C,Pm,Pt,Sw,Cw,K,PR,MCp,CCp,SGM)

if isempty(Pm)==0 && isempty(Pt)==0
    
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

Bw = SGM.Bwog(:,2);
Bo = SGM.Bwog(:,4);
Bg = SGM.Bwog(:,5);
Rs = SGM.Rs(:,2);

Bw1 = 0.5*(Bw(c) + Bw(r));
Bo1 = 0.5*(Bo(c) + Bo(r));
Bg1 = 0.5*(Bg(c) + Bg(r));
Rs1 = 0.5*(Rs(c) + Rs(r));
Cgsw1 = 0.5*(SGM.Cgsw(c) + SGM.Cgsw(r));
Cgso1 = 0.5*(SGM.Cgso(c) + SGM.Cgso(r));

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

Afw=Sat_cal(Swe(vP==1),1,1,as,aw)./Bw1(vP==1)./wmu1; 
Afp=Sat_cal(Swe(vP==1),1,1,as,aw)./Bw1(vP==1)./wmu1.*Cpe(vP==1); 
Afo=Sat_cal(Swe(vP==1),2,1,as,aw)./Bo1(vP==1)/mu(2); %oil
Afg=Sat_cal(Swe(vP==1),3,1,as,aw)./Bg1(vP==1)/mu(3);

Cfw=Sat_cal(Swe(vP==0),1,1,ts,tw)./Bw1(vP==0)./wmu0; %water
Cfp=Sat_cal(Swe(vP==0),1,1,ts,tw)./Bw1(vP==0)./wmu0.*Cpe(vP==0);
Cfo=Sat_cal(Swe(vP==0),2,1,ts,tw)./Bo1(vP==0)/mu(2); %oil
Cfg=Sat_cal(Swe(vP==0),3,1,ts,tw)./Bg1(vP==0)/mu(3);
%Cfp=Sat_tube(Swe(vP==0).*Cpe(vP==0),1,1,ts,tw)/mu(4); %polim

v1=find(vP==1);
v2=find(vP==0);
m1=size(vP,1);
vone=ones(m1,1);
% 
 %size([P])
% size([v1;v2])
% size(vone)
% size([Afw+Afo;Cfw+Cfo])

%ACL=sparse([v1;v2],vone,[Afg + (Rs1(vP==1)-Cgso1(vP==1)).*Afo - Cgsw1(vP==1).*Afw;Cfg + (Rs1(vP==0)-Cgso1(vP==0)).*Cfo - Cgsw1(vP==0).*Cfw],m1,1);
ACW=sparse([v1;v2],vone,[Afw;Cfw],m1,1);
ACO=sparse([v1;v2],vone,[Afo;Cfo],m1,1);
ACWCgs=sparse([v1;v2],vone,[-Cgsw1(vP==1).*Afw;-Cgsw1(vP==0).*Cfw],m1,1);
ACOCgs=sparse([v1;v2],vone,[(Rs1(vP==1)-Cgso1(vP==1)).*Afo;(Rs1(vP==0)-Cgso1(vP==0)).*Cfo],m1,1);
ACG=sparse([v1;v2],vone,[Afg;Cfg],m1,1);  %////////
ACP=sparse([v1;v2],vone,[Afp;Cfp],m1,1);

%ACL
a2cl=a2c.*(ACG + ACWCgs + ACOCgs);
a2cw=a2c.*ACW;
a2co=a2c.*ACO;
a2cwcgs=a2c.*ACWCgs;
a2cocgs=a2c.*ACOCgs;
a2cg=a2c.*ACG;
a2cp=a2c.*ACP;

A2CL=sparse(r,c,a2cl,n,m);
A2CW=sparse(r,c,a2cw,n,m);
A2CO=sparse(r,c,a2co,n,m);
A2CG=sparse(r,c,a2cg,n,m);
A2CP=sparse(r,c,a2cp,n,m);
A2CWCgs=sparse(r,c,a2cwcgs,n,m);
A2COCgs=sparse(r,c,a2cocgs,n,m);
else
    A2CL=A2C;
    A2CW=A2C;
    A2CO=A2C;
    A2CG=A2C;
    A2CP=A2C;
    A2CWCgs=A2C;
    A2COCgs=A2C;
end