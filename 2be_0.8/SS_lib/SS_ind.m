function [CR]=SS_ind(RC,na)

 v1=zeros(na,1);
 v1([RC.ACr;RC.AGr])=1;
 de=find(v1==1);
 
[rc1,rc2]=RC4IND(RC.Arc(:,1),RC.Arc(:,2));
 A_rc(:,1)=rc1;
 A_rc(:,2)=rc2;
 
 [C_r,C_c]=RC4IND(RC.Cr,RC.Cc);
 [G_r,G_c]=RC4IND(RC.Gr,RC.Gc);
 [AC_r,AC_c]=RC4IND(RC.ACr,RC.ACc);
 [AG_r,AG_c]=RC4IND(RC.AGr,RC.AGc);
 [CG_r,CG_c]=RC4IND(RC.CGr,RC.CGc);
 
 de1=de*2-1;
 de2=de*2;
 
 A_de=[de1,de2];

CR.A_rc=A_rc;
CR.C_r=C_r;
CR.C_c=C_c;
CR.G_r=G_r;
CR.G_c=G_c;

CR.AC_r=AC_r;
CR.AC_c=AC_c;
CR.AG_r=AG_r;
CR.AG_c=AG_c;

CR.CG_r=CG_r;
CR.CG_c=CG_c;

CR.A_de=A_de;

end

function [rc1,rc2]=RC4IND(r,c)
 r1=r*2-1;
 c1=c*2-1;
 
 r2=r*2-1;
 c2=c*2;
 
 r3=r*2;
 c3=c*2-1;
 
 r4=r*2;
 c4=c*2;
 
 rc1=[r1;r2;r3;r4];
 rc2=[c1;c2;c3;c4];
end