function dt=vibor_t2(dt,P,RC,dV,TL,W1,Won,Pw,na,PR,st,Ta,Sw,Sw2,dt0,Nl,ka1,va,vd,DL,W1D,WonD,nd)
Fc=PR.Fc;
Fc2=PR.Fc2;
Sc=PR.Sc;
Sc2=PR.Sc2;
% Sw=Sw([va,vd]);
% Sw2=Sw2([va,vd]);
dV1=dV(va);
dV2=dV(vd);

dP1=P(RC.Acr2(:,2))-P(RC.Acr2(:,1));
stnd=size(P,1)-nd;
dP2=P(stnd+RC.Dr2)+P(stnd+RC.Dc2);

PwNl=repmat(Pw,Nl,1);
PwNl=PwNl(ka1==1);
dPw1=P(Won)-PwNl;
dPw2=P(WonD(:,1))-PwNl(WonD(:,3));
v1=dP1>0;
TL=TL(RC.Acr2(:,1)+(RC.Acr2(:,2)-1)*na);

v2=dP2>0;
DL=DL(RC.Dr2+(RC.Dc2-1)*nd);

if dt==0
    if (min(Sw(va))>Sc) && (min(Sw(vd))>Sc2)
        MdS=abs((Sw2-Sw));
        nf=(isnan(MdS)==0);
        dS=max(MdS(nf==1));
      %  dS
        new_dt=0.005/dS;
        dt=new_dt*dt0;
      %  dt
    else
                
        dv=dV1(RC.Acr2(:,2)).*(v1==0)+dV1(RC.Acr2(:,1)).*v1;
        dt1=1./max(abs(TL.*dP1./dv));
        dt2=1./max(abs(W1.*dPw1./dV1(Won)));
        dt12=min([dt1,dt2])/Fc;
        
        dv=dV2(RC.Dr2).*(v2==0)+dV2(RC.Dc2).*v2;
        dt3=1./max(abs(DL.*dP2./dv));
        dt4=1./max(abs(W1D.*dPw2./dV2(WonD(:,1))));
        dt34=min([dt3,dt4])/Fc2;
        
        dt=min([dt12,dt34]);
    end;
    

    %kjhjkhj
end;

    if st+dt>Ta
        dt=Ta-st;
    end;