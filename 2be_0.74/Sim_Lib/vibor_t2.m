function dt=vibor_t2(dt,P,RC,dV,TL,W1,Won,Pw,na,PR,st,Ta,Sw,Sw2,dt0,Nl,WonM,va,vd,DL,W1D,WonD,nd,dV1,dV2)
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
dPw1=P(Won)-PwNl(WonM);
dPw2=P(WonD(:,1))-PwNl(WonD(:,3));
v1=dP1>0;
TL=TL(RC.Arc2(:,1)+(RC.Arc2(:,2)-1)*na);

v2=dP2>0;
Dl=zeros(0,1);
Dl=[Dl;DL(RC.Dr2+(RC.Dc2-1)*nd)];

    fl_swm=(min(Sw(va))>Sc);
if numel(vd)~=0
    fl_swd=(min(Sw(vd))>Sc2);
else
    fl_swd=-1;
end

if dt==0
    if fl_swm==1
        MdS=abs((Sw2(va)-Sw(va)));
        nf=(isnan(MdS)==0);
        dS=max(MdS(nf==1));
      %  dS
        new_dt=0.005/dS;
        dtm=new_dt*dt0;
      %  dt
    else
        dv=dV1(:,1).*(v1==0)+dV1(:,2).*v1;
        dt1=1./max(abs(TL.*dP1./dv));
        dt2=1./max(abs(W1.*dPw1./dVa(Won)));
        dt12=min([dt1,dt2])/Fc;
        dtm=min(dt12);
    end;
    
    if fl_swd==1
        MdS=abs((Sw2(vd)-Sw(vd)));
        nf=(isnan(MdS)==0);
        dS=max(MdS(nf==1));
      %  dS
        new_dt=0.005/dS;
        dtd=new_dt*dt0;
      %  dt
    else
        dv=dV2(:,1).*(v2==0)+dV2(:,2).*v2;
        dt3=1./max(abs(Dl.*dP2./dv));
        dt4=1./max(abs(W1D.*dPw2./dVd(WonD(:,1))));
        dt34=min([dt3;dt4])/Fc2;
        dtd=min(dt34);
    end;

   dt=min([dtm,dtd]);
end;

    if st+dt>Ta
        dt=Ta-st;
    end;