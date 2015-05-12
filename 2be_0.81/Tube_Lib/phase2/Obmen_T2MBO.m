function [A2CW,A2CO,A2CG,A2CORS,A2CP]=Obmen_T2MBO(A2C,Pm,Pt,K,PR,MCp,CCp,KWOG,CMP,va,vc)

if isempty(Pm)==0 && isempty(Pt)==0
    
mu=PR.mu;
mup=PR.mup;

[n,m]=size(A2C);
[r,c]=find(A2C);
A2C=A2C(:);
a2c=A2C(r+(c-1)*n).*K(r);

Bw1 = 0.5*(CMP.Bw(vc(c),2) + CMP.Bw(va(r),2));
Bo1 = 0.5*(CMP.Bo(vc(c),2) + CMP.Bo(va(r),2));
Bg1 = 0.5*(CMP.Bg(vc(c),2) + CMP.Bg(va(r),2));

vP=Pm(r,:)>=Pt(c,:);

Kwc=KWOG.w(vc(c));
Kwl=KWOG.w(va(r));
    
Koc=KWOG.o(vc(c));
Kol=KWOG.o(va(r));

Kgc=KWOG.g(vc(c));
Kgl=KWOG.g(va(r));

Cpc=MCp(r);
Cpl=CCp(c);

Kfw=Kwl.*vP(:,1)+Kwc.*(vP(:,1)==0);
Kfo=Kol.*vP(:,2)+Koc.*(vP(:,2)==0);
Kfg=Kol.*vP(:,3)+Koc.*(vP(:,3)==0);
Cpe=Cpc.*vP(:,1)+Cpl.*(vP(:,1)==0);

Rs = CMP.Rs(va(r),2).*vP(:,2)+CMP.Rs(vc(c),2).*(vP(:,2)==0);

wmu1=fusion(Cpe(vP(:,1)==1),mup);
wmu0=fusion(Cpe(vP(:,1)==0),mup);

Afw = Kfw(vP(:,1)==1)./Bw1(vP(:,1)==1)./wmu1; 
Afp = Kfw(vP(:,1)==1)./Bw1(vP(:,1)==1)./wmu1.*Cpe(vP(:,1)==1); 
Afo = Kfo(vP(:,2)==1)./Bo1(vP(:,2)==1)/mu(2);
Afg = Kfg(vP(:,3)==1)./Bg1(vP(:,3)==1)/mu(3);

Cfw = Kfw(vP(:,1)==0)./Bw1(vP(:,1)==0)./wmu0; 
Cfp = Kfw(vP(:,1)==0)./Bw1(vP(:,1)==0)./wmu0.*Cpe(vP(:,1)==0);
Cfo = Kfo(vP(:,2)==0)./Bo1(vP(:,2)==0)/mu(2);
Cfg = Kfg(vP(:,3)==0)./Bg1(vP(:,3)==0)/mu(3);

%Cfp=Sat_tube(Swe(vP==0).*Cpe(vP==0),1,1,ts,tw)/mu(4); %polim

v1w=find(vP(:,1)==1); v1o=find(vP(:,2)==1); v1g=find(vP(:,3)==1);
v2w=find(vP(:,1)==0); v2o=find(vP(:,2)==0); v2g=find(vP(:,3)==0); 
m1=size(vP,1);
vone=ones(m1,1);

ACW=sparse([v1w;v2w],vone,[Afw;Cfw],m1,1);
ACO=sparse([v1o;v2o],vone,[Afo;Cfo],m1,1);
ACORS=sparse([v1o;v2o],vone,[Rs(vP(:,2)==1).*Afo;Rs(vP(:,2)==0).*Cfo],m1,1);
ACG=sparse([v1g;v2g],vone,[Afg;Cfg],m1,1);  %////////
ACP=sparse([v1w;v2w],vone,[Afp;Cfp],m1,1);

a2cw=a2c.*ACW;
a2cors=a2c.*ACORS;
a2co=a2c.*ACO;
a2cg=a2c.*ACG;
a2cp=a2c.*ACP;

A2CW=sparse(r,c,a2cw,n,m);
A2CO=sparse(r,c,a2co,n,m);
A2CORS=sparse(r,c,a2cors,n,m);
A2CG=sparse(r,c,a2cg,n,m);
A2CP=sparse(r,c,a2cp,n,m);
else
    A2CW=A2C;
    A2CO=A2C;
    A2CORS=A2C;
    A2CG=A2C;
    A2CP=A2C;
end