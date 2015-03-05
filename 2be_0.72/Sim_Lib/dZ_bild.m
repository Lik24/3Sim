function [dZ]=dZ_bild(Z,RC,g,Ro)

dZ(1,:)=dz_nat(Z,RC.Arc(:,1),RC.Arc(:,2),g,Ro,RC.na);
dZ(2,:)=dz_nat(Z(RC.ACr),RC.Cr,RC.Cc,g,Ro,RC.nc);
dZ(3,:)=dz_nat(Z(RC.AGr),RC.Gr,RC.Gc,g,Ro,RC.ng);
dZ(4,:)=dz_nat(Z(RC.ADr),RC.Dr,RC.Dc,g,Ro,RC.nd);

end

function dZ=dz_nat(Z,r,c,g,Ro,n)
dZA=(Z(r).*Ro(1)-Z(c).*Ro(1))*g;    dZW=sparse(r,c,dZA,n,n);
dZA=(Z(r).*Ro(2)-Z(c).*Ro(2))*g;    dZO=sparse(r,c,dZA,n,n);
dZA=(Z(r).*Ro(3)-Z(c).*Ro(3))*g;    dZG=sparse(r,c,dZA,n,n);

dZ(1,1)={dZW};
dZ(1,2)={dZO};
dZ(1,3)={dZG};

end
