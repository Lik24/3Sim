function [A2CW,A2CO,A2CG,A2CWCgs,A2COCgs,A2CP]=Obmen_T2M(A2C,Pm,Pt,Sw,Cw,K,PR,MCp,CCp,SGM,CMP)

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

Bw1 = 0.5*(CMP.Bw(c,2) + CMP.Bw(r,2));
Bo1 = 0.5*(CMP.Bo(c,2) + CMP.Bo(r,2));
Bg1 = 0.5*(CMP.Bg(c,2) + CMP.Bg(r,2));
Rs1 = 0.5*(CMP.Rs(c,2) + CMP.Rs(r,2));
Cw1 = 0.5*(SGM.Cw(c) + SGM.Cw(r));
Cor1 = 0.5*(SGM.Cor(c) + SGM.Cor(r));
Co1 = 0.5*(SGM.Co(c) + SGM.Co(r));
Cg1 = 0.5*(SGM.Cg(c) + SGM.Cg(r));

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

Afw = Sat_cal(Swe(vP==1),1,1,as,aw)./Bw1(vP==1)./wmu1; 
Afp = Sat_cal(Swe(vP==1),1,1,as,aw)./Bw1(vP==1)./wmu1.*Cpe(vP==1); 
Afo = Sat_cal(Swe(vP==1),2,1,as,aw)./Bo1(vP==1)/mu(2); %oil
Afg = Cg1(vP==1).*Sat_cal(Swe(vP==1),3,1,as,aw)./Bg1(vP==1)/mu(3);

Cfw = Sat_cal(Swe(vP==0),1,1,ts,tw)./Bw1(vP==0)./wmu0; %water
Cfp = Sat_cal(Swe(vP==0),1,1,ts,tw)./Bw1(vP==0)./wmu0.*Cpe(vP==0);
Cfo = Sat_cal(Swe(vP==0),2,1,ts,tw)./Bo1(vP==0)/mu(2); %oil
Cfg = Cg1(vP==0).*Sat_cal(Swe(vP==0),3,1,ts,tw)./Bg1(vP==0)/mu(3);

%Cfp=Sat_tube(Swe(vP==0).*Cpe(vP==0),1,1,ts,tw)/mu(4); %polim

v1=find(vP==1);
v2=find(vP==0);
m1=size(vP,1);
vone=ones(m1,1);

ACW=sparse([v1;v2],vone,[SGM.Cg(v1).*Afw + SGM.Cpb(v1).*Rs1(v1).*Afo;SGM.Cg(v2).*Cfw + SGM.Cpb(v2).*Rs1(v2).*Cfo],m1,1);
ACO=sparse([v1;v2],vone,[Afo;Cfo],m1,1);
ACWCgs=sparse([v1;v2],vone,[Cw1(v1).*Afw;Cw1(v2).*Cfw],m1,1);
ACOCgs=sparse([v1;v2],vone,[(Cor1(v1).*Rs1(vP==1)+Co1(v1)).*Afo;(Cor1(vP==0).*Rs1(v2)+Co1(v2)).*Cfo],m1,1);
ACG=sparse([v1;v2],vone,[Afg;Cfg],m1,1);  %////////
ACP=sparse([v1;v2],vone,[Afp;Cfp],m1,1);

a2cw=a2c.*ACW;
a2co=a2c.*ACO;
a2cwcgs=a2c.*ACWCgs;
a2cocgs=a2c.*ACOCgs;
a2cg=a2c.*ACG;
a2cp=a2c.*ACP;

A2CW=sparse(r,c,a2cw,n,m);
A2CO=sparse(r,c,a2co,n,m);
A2CG=sparse(r,c,a2cg,n,m);
A2CP=sparse(r,c,a2cp,n,m);
A2CWCgs=sparse(r,c,a2cwcgs,n,m);
A2COCgs=sparse(r,c,a2cocgs,n,m);
else
    A2CW=A2C;
    A2CO=A2C;
    A2CG=A2C;
    A2CP=A2C;
    A2CWCgs=A2C;
    A2COCgs=A2C;
end