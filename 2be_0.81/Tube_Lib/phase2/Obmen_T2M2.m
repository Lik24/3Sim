function [A2CW,A2CO,A2CP]=Obmen_T2M2(A2C,Pm,Pt,K,PR,MCp,CCp,KWOG,CMP,va,vc)

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

Bw1 = 0.5*(CMP.Bw(vc(c),2) + CMP.Bw(va(r),2));
Bo1 = 0.5*(CMP.Bo(vc(c),2) + CMP.Bo(va(r),2));

vP=Pm(r,:)>=Pt(c,:);

% vP
% [Pm(r),Pt(c)]

 Kwc=KWOG.w(vc(c));
 Kwl=KWOG.w(va(r));
    
 Koc=KWOG.o(vc(c));
 Kol=KWOG.o(va(r));
      
 Kfw=Kwl.*vP(:,1)+Kwc.*(vP(:,1)==0);
 Kfo=Kol.*vP(:,2)+Koc.*(vP(:,2)==0);

Cpl=MCp(r);
Cpc=CCp(c);

Cpe=Cpl.*vP(:,1)+Cpc.*(vP(:,1)==0);

wmu1=fusion(Cpe(vP(:,1)==1),mup);
wmu0=fusion(Cpe(vP(:,1)==0),mup);

Afw = Kfw(vP(:,1)==1)./Bw1(vP(:,1)==1)./wmu1; 
Afp = Kfw(vP(:,1)==1)./Bw1(vP(:,1)==1)./wmu1.*Cpe(vP(:,1)==1); 
Afo = Kfo(vP(:,2)==1)./Bo1(vP(:,2)==1)/mu(2); %oil

Cfw = Kfw(vP(:,1)==0)./Bw1(vP(:,1)==0)./wmu0; %water
Cfp = Kfw(vP(:,1)==0)./Bw1(vP(:,1)==0)./wmu0.*Cpe(vP(:,1)==0);
Cfo = Kfo(vP(:,2)==0)./Bo1(vP(:,2)==0)/mu(2); %oil

v1w=find(vP(:,1)==1); v1o=find(vP(:,2)==1);
v2w=find(vP(:,1)==0); v2o=find(vP(:,2)==0);
m1=size(vP,1);
vone=ones(m1,1);

ACW=sparse([v1w;v2w],vone,[Afw;Cfw],m1,1);
ACO=sparse([v1o;v2o],vone,[Afo;Cfo],m1,1);
ACP=sparse([v1w;v2w],vone,[Afp;Cfp],m1,1);

a2cw=a2c.*ACW;
a2co=a2c.*ACO;
a2cp=a2c.*ACP;

A2CW=sparse(r,c,a2cw,n,m);
A2CO=sparse(r,c,a2co,n,m);
A2CP=sparse(r,c,a2cp,n,m);

else
    A2CW=A2C;
    A2CO=A2C;
    A2CP=A2C;
end