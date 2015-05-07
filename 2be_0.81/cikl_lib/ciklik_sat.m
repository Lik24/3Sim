function Qc=ciklik_sat(Sw_old,PR,qe,rc)

dxdy=1;
Kwi=Sat_cal(Sw_old,1,1,PR.as,PR.aw)/PR.mu(1);
Koi=Sat_cal(Sw_old,2,1,PR.as,PR.aw)/PR.mu(2);

f=Kwi./(Kwi+Koi);
df=f(rc(:,1))-f(rc(:,2));
df=sparse(rc(:,1),rc(:,2),df);

Qc=PR.BetCikl*dxdy*(df*qe);