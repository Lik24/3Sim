function CrDATA=CrackProp(DATA,PR,NT)

H=DATA.gH;
dC=PR.dl;

kc=PR.Kc;
dh=PR.dh;

for l=1:PR.Nl
  nt=NT{l};
  KC(l)={kc*ones(size(nt))};
  DH(l)={dh*ones(size(nt))};
end

CrDATA.H=H;
CrDATA.dC=dC;
CrDATA.KC=KC;
CrDATA.DH=DH;