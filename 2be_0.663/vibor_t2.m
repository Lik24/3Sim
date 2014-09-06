function dt=vibor_t2(dt,P,RC,dV,TL,W1,Won,Pw,na,PR,st,Ta)
Fc=PR.Fc;

dP=P(RC.Acr2(:,2))-P(RC.Acr2(:,1));
dPw=P(Won)-Pw;
v1=dP>0;
TL=TL(RC.Acr2(:,1)+(RC.Acr2(:,2)-1)*na);

if dt==0
   dv=dV(RC.Acr2(:,2)).*(v1==0)+dV(RC.Acr2(:,1)).*v1;
   dt1=1./max(abs(TL.*dP./dv));
%    TL.*dP
%    dv
   
  dt2=1./max(abs(W1.*dPw./dV(Won)));
%    dt1
%    dt2
%   Fc
  dt=min([dt1,dt2])/Fc;

  if st+dt>Ta
      dt=Ta-st;
  end;
  %kjhjkhj
end;