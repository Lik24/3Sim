function [A2CL,A2CW,A2CP,A2CG]=Obmen_T2M_3(a2c,Pm,Pt,Sw,Cw,ts,tw,as,aw,mu,MCp,CCp,r,c,kfw,kfo)

if isempty(r)==0
    m=size(Cw,1);
    
    kaf=Sat_cal(Sw,1,1,as,aw); %water
    kaof=Sat_cal(Sw,2,1,as,aw); %oil
    
    kcf=Sat_tube(Cw,1,1,ts,tw); %water
    kcof=Sat_tube(Cw,2,1,ts,tw); %oil
    
    vP=Pm(r)>=Pt(c);
    
    Kafc=kaf(r);
    Kafl=kcf(c);
    
    Kaofc=kaof(r);
    Kaofl=kcof(c);
    
    Kwe=Kafc.*vP+Kafl.*(vP==0);
    Koe=Kaofc.*vP+Kaofl.*(vP==0);
    
    Cpc=MCp(r);
    Cpl=CCp(c);
    
    Cpe=Cpc.*vP+Cpl.*(vP==0);
    
    Afo=Koe(vP==1)./mu(2);
    Cfo=Koe(vP==0)./mu(2);
    
    Afw=Kwe(vP==1).*((1-Cpe(vP==1))./mu(1)+Cpe(vP==1)/mu(4));
    Cfw=Kwe(vP==0).*((1-Cpe(vP==0))./mu(1)+Cpe(vP==0)/mu(4));
    
    Afp=Kwe(vP==1).*Cpe(vP==1)/mu(4);
    Cfp=Kwe(vP==0).*Cpe(vP==0)/mu(4);
    
    v1=find(vP==1);
    v2=find(vP==0);
    vone=ones(m,1);
    
    ACO=sparse([v1;v2],vone,[Afo;Cfo],m,1);
    ACW=sparse([v1;v2],vone,[Afw;Cfw],m,1);
    ACP=sparse([v1;v2],vone,[Afp;Cfp],m,1);
    
    A2CL=a2c.*(ACW+ACO);
    A2CW=a2c.*ACW;
    A2CP=a2c.*ACP;
    
    A2CG=1;
else
    
    A2CL=a2c;
    A2CW=a2c;
    A2CP=a2c;
    
    A2CG=1;
end;