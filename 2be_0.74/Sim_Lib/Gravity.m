function [Gr,Grw]=Gravity(TL,TW,CL,CW,GL,GW,DL,DW,TG,A,dZ)

TO=TL-TW;
CO=CL-CW;
GO=GL-GW;
DO=DL-DW;
[grw_A,gro_A]=degrav(dZ(1,:),TW,TO);
[grw_C,gro_C]=degrav(dZ(2,:),CW',CO');
[grw_G,gro_G]=degrav(dZ(3,:),GW,GO);
[grw_D,gro_D]=degrav(dZ(4,:),DW,DO);

Gr=[gro_A;gro_C;gro_G;gro_D]+A.*[grw_A;grw_C;grw_G;grw_D];
Grw=[grw_A;grw_C;grw_G;grw_D];
end

function [Grw,Gro]=degrav(dZ,TW,TO)
dZW=dZ{1};
dZO=dZ{2};
dZG=dZ{3};

Grw=sum(dZW.*TW,2);
Gro=sum(dZO.*TO,2);
%Grg=sum(dZG.*TG,2);
end
