function T2=Termal(T1,P,CSw,Mp,Mc,Mg,NTG,LM,LC1,LG1,LA2C,LA2G,PR,RC,dVCG,dt,Won,WonM,WonV,WonG,TeW,Qm,Qc,Qg,...
    C2G,TL,CL,GL,TW,CW,GW,A2C,A2CW,A2G,A2GW,Wf,ndt,t,T0,S,BZ)

La=PR.lam;
Cp=PR.Cp;
Ro=PR.Ro;

na=RC.na;
nc=RC.nc;
ng=RC.ng;

va=1:na;
vc=na+1:na+nc;
vg=na+nc+1:na+nc+ng;

Sw=CSw(va);  Cw=CSw(vc);   Gw=CSw(vg);
Pm=P(va);    Pc=P(vc);     Pg=P(vg);

bm=(Cp(1).*Ro(1).*Sw.*Mp+Cp(2).*Ro(2).*(1-Sw).*Mp+Cp(4).*Ro(4).*(1-Mp)).*NTG+Cp(5).*Ro(5).*(1-NTG);
bm=bm./Mp;
bc=Cp(1).*Ro(1).*Cw.*Mc+Cp(2).*Ro(2).*(1-Cw).*Mc+Cp(4).*Ro(4).*(1-Mc);
bg=Cp(1).*Ro(1).*Gw.*Mg+Cp(2).*Ro(2).*(1-Gw).*Mg+Cp(4).*Ro(4).*(1-Mg);
%ndt
ndt=1;
%Cp(:)=0;
[LA,LA0,LA_CND]=Potok_tepl(LM,Sw,La,NTG,Mp,RC.Arc(:,1),RC.Arc(:,2),TL,TW,Pm,Cp,Ro,ndt,dt);
[LC,LC0,LC_CND]=Potok_tepl(LC1,Cw,La,1,Mc,RC.Cr,RC.Cc,CL,CW,Pc,Cp,Ro,ndt,dt);
[LG,LG0,LG_CND]=Potok_tepl(LG1,Gw,La,1,Mg,RC.Gr,RC.Gc,GL,GW,Pg,Cp,Ro,ndt,dt); % ��������� ����������


%  full(LC')
%  fghj
[A2CLi,A2CLo]=Obmen_Termo(LA2C,Sw,Cw,La,NTG,Mp,Mc,RC.Cr,RC.Cc,na,nc,Pm,Pc,A2C,A2CW,Cp,Ro,ndt,dt);
[A2GLi,A2GLo]=Obmen_Termo(LA2G,Sw,Gw,La,NTG,Mp,Mg,RC.Gr,RC.Gc,na,ng,Pm,Pg,A2G,A2GW,Cp,Ro,ndt,dt);

%Qm

[WAcnvi,WAcnvI,WAcnvo,WAcnd,WAcnd2]=Well_Term(Qm,TeW(WonM),Won,Ro,Cp,na,ndt,La,Wf,dt);
[WCcnvi,WCcnvI,WCcnvo,WCcnd,WCcnd2]=Well_Term(Qc,TeW(WonV(:,3)),WonV(:,1),Ro,Cp,nc,ndt,La,WonV(:,2),dt);
[WGcnvi,WGcnvI,WGcnvo,WGcnd,WGcnd2]=Well_Term(Qg,TeW(WonG(:,3)),WonG(:,2),Ro,Cp,ng,ndt,La,WonG(:,1),dt);

[bmp]=MPSS(Cp.*Ro,La,t*dt,T1,T0,S,BZ);

VB=[bm;bc;bg].*dVCG;

    for i=1:ndt
% LA(1:1000)
% LA_CND(1:1000)
% fgh
        A1=(LA+LA_CND+sparse(1:na,1:na,-VB(va)+WAcnvI+WAcnvo-WAcnd+sum(A2CLi,2)+sum(A2CLo,2),na,na));
        C1=(LC+LC_CND+sparse(1:nc,1:nc,-VB(vc)+WCcnvI+WCcnvo-WCcnd-sum(A2CLo,1)'-sum(A2CLi,1)',nc,nc));
        G1=(LG+LG_CND+sparse(1:ng,1:ng,-VB(vg)+WGcnvI+WGcnvo-WGcnd,ng,ng));
        
        A=[A1,-A2CLo,-A2GLo;A2CLi',C1,C2G;A2GLi',C2G',G1];

        T1(:,1)=(-VB.*T1-[WAcnvi;WCcnvi;WGcnvi]-[WAcnvo;WCcnvo;WGcnvo].*T1-[LA0,LC0,LG0]'.*T1...
            -[sum(A2CLo,2)-sum(A2GLi,2);-sum(A2CLi,1)';-sum(A2GLi,1)'].*T1...
            -[WAcnd2;WCcnd2;WGcnd2]+[bmp;zeros(nc+ng,1)])'/A;
       % T1(:,1)=T1.*(T1>=40).*(T1<=100)+40.*(T1<40)+100.*(T1>100);
%         T1
%         ljkbh
    end;

T2=T1;
% sum(VB(vc).*(T2(vc)-T11(vc)))