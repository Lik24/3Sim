function [ndt,j_ndt,fl1]=vibor_t(ndt,fl,Pc,Pg,Pw,PR,RC,dt,Sw1,Sw2,CL,GL,dVc,dVg,W1C,W1G,WoC,WoG,i,j_ndt)

Fc=PR.Fc2;
Sc=PR.Sc2;

dPc=Pc(RC.Cc2)-Pc(RC.Cr2);
dPg=Pg(RC.Gc2)-Pg(RC.Gr2);
dPwc=Pc(WoC(:,1))-Pw(WoC(:,3));
dPwg=Pg(WoG(:,1))-Pw(WoG(:,3));

v1=dPc>0;
v2=dPg>0;

if  fl==0
  MdS=abs(Sw2-Sw1);
  Sw_ch=Sw1(MdS~=0);
  if isempty(Sw_ch)==1
      Sw_ch=1;
  end
  %i=1;
  if (min(Sw_ch)>=Sc) && (i~=1)
    MdS=abs(Sw2-Sw1);
    nf=(isnan(MdS)==0);
    dS=max(MdS(nf==1));
    new_dt=dS/0.01;
    ndt=ceil(ndt*new_dt);    
  else
    dvc=dVc(RC.Cc2).*(v1==0)+dVc(RC.Cr2).*v1;
    dvg=dVg(RC.Gc2).*(v2==0)+dVg(RC.Gr2).*v2;
    
    dt1=1./max(abs([CL.*dPc./dvc;GL.*dPg./dvg]));
    ndt1=dt/dt1;    
    
    if isempty(WoC(:,1))~=1 || isempty(WoG(:,1))~=1 
     WdS=abs(Sw2(WoC(:,1))-Sw1(WoC(:,1)));
     QC=W1C.*dPwc./dVc(WoC(:,1));
     QC=QC(WdS~=0);
     dt2=1./max(abs([QC;W1G.*dPwg./dVg(WoG(:,1))]));
     ndt2=dt/dt2;
    else 
     ndt2=0;   
    end;
%     ndt1
%     ndt2
%     sdf
    ndt=ceil(max([ndt1;ndt2])*Fc);
  end;
else
   ndt=1;
end;

j_old=j_ndt;
j_ndt=j_ndt+1/ndt;

fl1=j_ndt<1;

if fl1==0
    ndt=1/(1-j_old);
end;