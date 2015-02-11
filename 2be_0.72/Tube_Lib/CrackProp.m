function CrDATA=CrackProp(DATA,PR,NT)

H=DATA.gH;
dC=PR.dl;

kc=PR.Kc;
dh=PR.dh;
kd=PR.Kd;
kg=PR.Kg;
hg=PR.Hg;

alp_C=PR.Alp_C;
KC={kc};
for i=1:size(NT,2)%PR.Nl
  %nt=NT{l};
  DH(i)={dh*ones(size(NT))};
  GH(i)={hg*ones(size(NT))};
  KC(i)={kc};
end

for l=1:PR.Nl
  %nt=NT{l};
  DH(l)={dh*ones(size(NT))};
  GH(l)={hg*ones(size(NT))};

  KG(l)={kg*ones(size(NT))};
  KD(l)={kd*ones(size(NT))};
end

CrDATA.H=H;
CrDATA.dC=dC;

CrDATA.DH=DH;
CrDATA.GH=GH;

CrDATA.KC=KC;
CrDATA.alp_C=alp_C;

CrDATA.KD=KD;
CrDATA.KG=KG;
CrDATA.Xk=kd;