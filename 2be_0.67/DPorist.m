function nd=DPorist(XY,Nl)
nd=cell(Nl,1);
n=size(XY,1);

for l=1:Nl
    LN=n*(l-1);
  for i=1:1
    nd_sl={(1:n)+LN};
  end
  nd(l)={nd_sl};  
end