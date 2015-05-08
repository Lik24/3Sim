function dt=vibor_t2(dt,P,RC,dV,TL,W1,WELL,Pw,na,PR,st,Ta,Sw,Sw2,dt0,Nl,va,vd,DL,W1D,nd,dV1,dV2)
Fc=PR.Fc;
Fc2=PR.Fc2;
Sc=PR.Sc;
Sc2=PR.Sc2;
% Sw=Sw([va,vd]);
% Sw2=Sw2([va,vd]);
dVa=dV(va);
dVd=dV(vd);

dP1=P(RC.Arc2(:,2))-P(RC.Arc2(:,1));
stnd=size(P,1)-nd;
dP2=P(stnd+RC.Dc2,1)-P(stnd+RC.Dr2,1);

PwNl=repmat(Pw,Nl,1);
%PwNl=PwNl(ka1==1);
dPw1=P(WELL.Won(:,1))-PwNl(WELL.Won(:,3));
dPw2=P(WELL.WonD(:,1))-PwNl(WELL.WonD(:,3));
v1=dP1>0;
TL=TL(RC.Arc2(:,1)+(RC.Arc2(:,2)-1)*na);

v2=dP2>0;
Dl=zeros(0,1);
Dl=[Dl;DL(RC.Dr2+(RC.Dc2-1)*nd)];

if numel(vd)~=0
    fl_sw=(min(Sw(va))>Sc)*(min(Sw(vd))>Sc2);
else
    fl_sw=(min(Sw(va))>Sc);
end

if dt==0
    if fl_sw==1
        MdS=abs((Sw2-Sw));
        nf=(isnan(MdS)==0);
        dS=max(MdS(nf==1));
      %  dS
        new_dt=0.005/dS;
        dt=new_dt*dt0;
      %  dt
    else
                
        dv=dV1(:,1).*(v1==0)+dV1(:,2).*v1;
        dt1=1./max(abs(TL.*dP1./dv));
        dt2=1./max(abs(W1.*dPw1./dVa(WELL.Won(:,1))));
        dt12=min([dt1,dt2])/Fc;
        
        dv=dV2(:,1).*(v2==0)+dV2(:,2).*v2;
        dt3=1./max(abs(Dl.*dP2./dv));
        dt4=1./max(abs(W1D.*dPw2./dVd(WELL.WonD(:,1))));
        dt34=min([dt3,dt4])/Fc2;
        
        dt=min([dt12;dt34(:,1)]);
    end;
    

    %kjhjkhj
end;

    if st+dt>Ta
        dt=Ta-st;
    end;