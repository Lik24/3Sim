function CrDATA=CrackProp(DATA,PR,NT)

H=DATA.gH;
dC=PR.dl;

kc=PR.Kc;
dh=PR.dh;
kd=PR.Kd;
kg=PR.Kg;
hg=PR.Hg;

for l=1:PR.Nl
  nt=NT{l};
  DH(l)={dh*ones(size(nt))};
  GH(l)={hg*ones(size(nt))};

  KC(l)={kc*ones(size(nt))};
  KG(l)={kg*ones(size(nt))};
  KD(l)={kd*ones(size(nt))};
end

CrDATA.H=H;
CrDATA.dC=dC;

CrDATA.DH=DH;
CrDATA.GH=GH;

CrDATA.KC=KC;
CrDATA.KD=KD;
CrDATA.KG=KG;
CrDATA.Xk=kd;