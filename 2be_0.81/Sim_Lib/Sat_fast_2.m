function [Sw,Cp,ndt,Q1,Q2,Qm,Qd,tmp,ndti]=Sat_fast_2(Sw,Cp,RC,TC,TG,A2C,A2G,A2D,D2C,D2G,Pi,PR,...
    ndt0,Won,Uf,dt,dV,Pw,WonG,CpW,WonC,Nl,CR_cr,Qz,Qf,Pi0,TL,W1,TW,W6,TP,...
    W7,L,Lc,Lg,Ke,Cws,Cwp,BLGY_GIM,Qzm1,WonM,nw,b1gm,b1gd,GYData,...
    Clp,ka1,W1D,W6D,W7D,A2BW,A2BP,D2BW,D2BP,DL,DW,DP,WoD,A2DW,A2DP,A2DL,BB,A2BL,D2BL,b1gb,Grw,dZ,Mp,Bwo,P0)

na=RC.na;
nc=RC.nc;
ng=RC.ng;
nd=RC.nd;
nb=RC.nb;
vad=RC.ADr;
b1wb=sparse(nb,1);

va=1:na;
vc=na+1:na+nc;
vg=na+nc+1:na+nc+ng;
vd=na+nc+ng+1:na+nc+ng+nd;
vb=na+nc+ng+nd+1:na+nc+ng+nd+nb;

PwNl=repmat(Pw,Nl,1);
%PwNl=PwNl(ka1==1);

     v1a=zeros(na,1);
     v1a(RC.ACr)=1;
     r1a=find(v1a==1);
     
     v2a=zeros(na,1);
     v2a(RC.AGr)=1;
     r2a=find(v2a==1);
 
     v1d=zeros(nd,1);
     v1d(RC.DCr)=1;
     r1d=find(v1d==1);
     
     v2d=zeros(nd,1);
     v2d(RC.DGr)=1;
     r2d=find(v2d==1);
% aw1=sum(SCw(vc).*dV(vc));

 if isempty(RC.Cr)==0 || isempty(RC.Gr)==0
     [Bwc,Bwg,Blc,Blg,Bwcd,Bwgd,Blcd,Blgd,Sw(vc),Sw(vg),ndt,Q1,Q2,Qm,Qd,~,ndti]=fun1(RC,Pi,Sw,Cp,PR,TC,TG,A2C,...
         A2G,A2D,D2C,D2G,WonC,WonG,Uf,CpW,Pw,dt,dV,CR_cr,Qz,Qf,ndt0,Pi0,L,Lc,Lg,Ke,dZ,Mp,Bwo,P0);
     %  [Bc,Bg,SCw(vc),SCw(vg),ndt,Q1,Q2,Qm,dSS]=fun2(RC,Pi,SCw,SCp,PR,TC,TG,A2C,...
     %     A2G,WonC,WonG,Uf,CpW,Pw,dt,dV,CR_cr,Qz,Qf,ndt0,Pi0,L,Lc,Lg,Ke,CR,TM);
     
     b1wm=sparse(Won(:,1),ones(1,size(Won,1)),-W1.*Pw(Won(:,3)),na,1);
     b1wd=sparse(WoD(:,1),ones(1,size(WoD,1)),-W1D.*Pw(WoD(:,3)),nd,1);
     
     W2M=sparse(WonM,Won(:,1),W1,nw,na);
     W2D=sparse(WoD(:,3),WoD(:,1),W1D,nw,nd);
     W2B=sparse(nw,nb);

    Qf=Qzm1<0;

     b1wm=b1wm.*(sum(W2M(Qf==0,:),1)~=0)';
     b1wd=b1wd.*(sum(W2D(Qf==0,:),1)~=0)';
     
     WM1=[W2M,W2D,W2B];
     WM2=WM1';
     W3vec=sparse(WonM,1,W1,nw,1)+sparse(WoD(:,3),1,W1D,nw,1);
     WM3=-sparse(1:nw,1:nw,W3vec,nw,nw);
  
     bl=[b1wm;b1wd;b1wb]+BLGY_GIM;
     Blc1=zeros(na,1);      Blc1(RC.ACr)=Blc;
     Bwc1=zeros(na,1);      Bwc1(RC.ACr)=Bwc;
     
     Blcd1=zeros(nd,1);     Blcd1(RC.DCr)=Blcd;
     Bwcd1=zeros(nd,1);     Bwcd1(RC.DCr)=Bwcd;
     
     blm=bl(va)-Blc1/dt-sparse(r2a,ones(sum(v2a),1),Blg,na,1)/dt;%sparse(r1a,ones(sum(v1a),1),Blc(),na,1)
     bld=bl(vd-nc-ng)-Blcd1/dt-sparse(r2d,ones(sum(v2d),1),Blgd,nd,1)/dt;
     Bl=[blm;bld;bl(vb-na-nd)];
     
     A1=TL-sparse(Won(:,1),Won(:,1),W1,na,na)-sparse(1:na,1:na,Clp(va)+sum(b1gm(:,1:2),2)+sum(A2DL,2),na,na);
     D1=DL-sparse(WoD(:,1),WoD(:,1),W1D,nd,nd)-sparse(1:nd,1:nd,Clp(vd)+sum(b1gd(:,1:2),2)+sum(A2DL,1)',nd,nd);                       %������� ����. ��� ������� ���.
     B1=BB-sparse(1:nb,1:nb,sum(A2BL,1)+sum(D2BL,1)+Clp(vb)'+b1gb',nb,nb);
        
     AM=[A1,   A2DL, A2BL;
        A2DL',  D1,  D2BL;
        A2BL', D2BL', B1];
     
     WM1=WM1(Qf~=0,:);
     WM2=WM2(:,Qf~=0);
     WM3=WM3(Qf~=0,Qf~=0);
     
     Qzm=zeros(size(Qzm1));
     Qzd=zeros(size(Qzm1));
     Qzm1(CR_cr(1,1).won(:,3))=0;
     Qzm1(CR_cr(1,2).won(:,3))=0;
     
     Qzm(CR_cr(1,1).won(:,3))=(Qm(CR_cr(1,1).won(:,3),1)+Qm(CR_cr(1,1).won(:,3),2))/dt;
     Qzd(CR_cr(1,2).won(:,3))=(Qd(CR_cr(1,2).won(:,3),1)+Qd(CR_cr(1,2).won(:,3),2))/dt;
     Qzm1=Qzm1+Qzm+Qzd;
     
     Pt=[Bl',Qzm1(Qf~=0)']/[AM,WM2;WM1,WM3];
     
     Pi(va)=Pt(va)';
     Pi(vd)=Pt(vd-na-nc-ng)';
     Pi(vb)=Pt(vb-na-nc-ng-nd)';
     PwNl(Qf~=0)=Pt(na+nd+nb+1:end);
     
     %��������
     Bpc=zeros(size(Bwc));
     Bpg=zeros(size(Bwg));
     Bpcd=zeros(size(Bwcd));
     Bpgd=zeros(size(Bwgd));

     Qm=QBild(W1,W6,W7,Pi(va,1),Uf(WonM),Won(:,1),dt,PwNl(WonM),WonM,nw);%PwNl
     Qd=QBild(W1D,W6D,W7D,Pi(vd,1),Uf(WoD(:,3)),WoD(:,1),dt,PwNl(WoD(:,3)),WoD(:,3),nw);
 else
     Bwc=zeros(nc,1);
     Bwg=zeros(ng,1);
     
     Bpc=zeros(size(Bwc));
     Bpg=zeros(size(Bwg));
    
     Bwcd=zeros(nc,1);
     Bwgd=zeros(ng,1);
     
     Bpcd=zeros(size(Bwcd));
     Bpgd=zeros(size(Bwgd));
     
     ndt=1;
     ndti=1;
     Q1=zeros(nw,5);
     Q2=zeros(nw,5);
     Qm=zeros(0,5);
     Qd=zeros(0,5);
     
     Blc1=zeros(na,1);     
     Bwc1=zeros(na,1);     
     
     Blcd1=zeros(nd,1);    
     Bwcd1=zeros(nd,1);   
 end;

b_A2B=Soed2B(A2BW,A2BP,Pi,na,nb,va,vb);     % ����� ��� � ��������� ��������
b_D2B=Soed2B(D2BW,D2BP,Pi,nd,nb,vd,vb);     % ����� ������ � ��������� ��������

     bwm=sparse(Won(:,1),ones(1,size(Won,1)),-W6.*(Pi(Won(:,1))-PwNl(WonM)),na,1)...
         -b1gm(:,3).*(Pi(va)-GYData.GY_Pxy)-b1gm(:,4).*(Pi(va)-GYData.GY_Pz)-b_A2B(:,1);
     
     bwd=sparse(WoD(:,1),ones(1,size(WoD,1)),-W6D.*(Pi(vd(WoD(:,1)))-PwNl(WoD(:,3))),nd,1)...
         -b1gd(:,3).*(Pi(vd)-GYData.GY_Pxy(vad))-b1gd(:,4).*(Pi(vd)-GYData.GY_Pz(vad))-b_D2B(:,1);
    
     bpm=sparse(Won(:,1),ones(1,size(Won,1)),-W7.*(Pi(Won(:,1))-PwNl(WonM)),na,1)...
         -b1gm(:,5).*(Pi(va)-GYData.GY_Pxy)-b1gm(:,6).*(Pi(va)-GYData.GY_Pz)-b_A2B(:,2);
     
     bpd=sparse(WoD(:,1),ones(1,size(WoD,1)),-W7D.*(Pi(vd(WoD(:,1)))-PwNl(WoD(:,3))),nd,1)...
         -b1gd(:,5).*(Pi(vd)-GYData.GY_Pxy(vad))-b1gd(:,6).*(Pi(vd)-GYData.GY_Pz(vad))-b_D2B(:,2);
     
     bw=[bwm;bwd];
     bp=[bpm;bpd];
%sparse([r1a;r1d+na],ones(sum(v1a)+sum(v1d),1),[Bwc;Bwcd],na+nd,1)/dt
     Bw=bw+[Bwc1;Bwcd1]/dt+sparse([r2a;r2d+na],ones(sum(v2a)+sum(v2d),1),[Bwg;Bwgd],na+nd,1)/dt...
         -Cwp([va,vd]).*(Pi([va,vd])-Pi0([va,vd]))+Grw([va,vd]);
     Bp=bp+sparse([r1a;r1d+na],ones(sum(v1a)+sum(v1d),1),[Bpc;Bpcd],na+nd,1)/dt+sparse([r2a;r2d+na],ones(sum(v2a)+sum(v2d),1),[Bpg;Bpgd],na+nd,1)/dt...
         -Cp([va,vd]).*Cwp([va,vd]).*(Pi([va,vd])-Pi0([va,vd]))+Cp([va,vd]).*Grw([va,vd]);
     tmp=sum(Bw);    

     AM2=[TW-sparse(1:na,1:na,sum(A2DW,2)),A2DW;
          A2DW',DW-sparse(1:nd,1:nd,sum(A2DW,1))];
     AM3=[TP,A2DP;A2DP',DP];
     
     Sw_old=Sw([va,vd]);
     Sw([va,vd])=Sw([va,vd])+dt*(AM2*Pi([va,vd])+Bw)./Cws([va,vd]);

     Cp([va,vd])=Sw_old.*Cp([va,vd])+dt*(AM3*Pi([va,vd])+Bp)./Cws([va,vd]);
     
     vad=[va,vd];
     v0=Sw(vad)==0;
     Cp(vad(v0==0))=Cp(vad(v0==0))./Sw(vad(v0==0));

    %Sw=Sw.*(Sw>=PR.aw(4)).*(Sw<=1-PR.aw(5))+(1-PR.aw(5)).*(Sw>1-PR.aw(5))+PR.aw(4).*(Sw<PR.aw(4));
    Swr=PR.Swr;
    Sor=PR.Sor;
     
    Sw=Sw.*(Sw>=Swr).*(Sw<=(1-Sor))+(1-Sor).*(Sw>(1-Sor))+Swr.*(Sw<Swr);
    %Sw=Sw.*(Sw>=0).*(Sw<=1)+(Sw>1); 
    Cp=Cp.*(Cp>=0).*(Cp<=1)+(Cp>1);
end

function b_A2B=Soed2B(A2BW,A2BP,Pi,na,nb,va,vb)
[r1,c1]=find(A2BW);
dP=(Pi(va(r1))-Pi(vb(c1)));
A2B=A2BW(r1+(c1-1)*na).*dP;
A2B1=sparse(r1,c1,A2B,na,nb);
b_A2B(:,1)=sum(A2B1,2);

A2B=A2BP(r1+(c1-1)*na).*dP;
A2B1=sparse(r1,c1,A2B,na,nb);
b_A2B(:,2)=sum(A2B1,2);
end