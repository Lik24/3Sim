function dt=vibor_t2(dt,P,RC,dV,TL,W1,Won,Pw,na,PR,st,Ta,Sw,Sw2,dt0,Nl,ka1,DL,W1D,WonD,nd)
Fc=PR.Fc;
Sc=PR.Sc2;

dP=P(RC.Acr2(:,2))-P(RC.Acr2(:,1));

PwNl=repmat(Pw,Nl,1);
PwNl=PwNl(ka1==1);
dPw=P(Won)-PwNl;
v1=dP>0;
TL=TL(RC.Acr2(:,1)+(RC.Acr2(:,2)-1)*na);


if dt==0
    if (min(Sw)>Sc)
        MdS=abs((Sw2-Sw));
        nf=(isnan(MdS)==0);
        dS=max(MdS(nf==1));
      %  dS
        new_dt=0.001/dS;
        dt=new_dt*dt0;
      %  dt
    else
        
        
        dv=dV(RC.Acr2(:,2)).*(v1==0)+dV(RC.Acr2(:,1)).*v1;
        dt1=1./max(abs(TL.*dP./dv));
        dt2=1./max(abs(W1.*dPw./dV(Won)));
        
        dt=min([dt1,dt2])/Fc;
    end;
    

    %kjhjkhj
end;

    if st+dt>Ta
        dt=Ta-st;
    end;