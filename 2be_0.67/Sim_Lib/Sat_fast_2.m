function [Sw,Cp,ndt,Q1,Q2,Qm,tmp]=Sat_fast_2(Sw,Cp,RC,TC,TG,A2C,A2G,Pi,PR,...
    ndt0,Won,Uf,dt,dV,Pw,WonG,CpW,WonC,Nl,CR_cr,Qz,Qf,Pi0,TL,W1,TW,W6,TP,...
    W7,L,Lc,Lg,Ke,Cws,Cwp,BLGY_GIM,Qzm1,WonM,nw,b1gm,b1gd,GYData,...
    Clp,ka1,W1D,W6D,W7D,A2BW,A2BP,D2BW,D2BP,DL,DW,DP,WoD,A2DW,A2DP)

na=RC.na;
nc=RC.nc;
ng=RC.ng;
nd=RC.nd;
nb=RC.nb;
vad=RC.ADr;

va=1:na;
vc=na+1:na+nc;
vg=na+nc+1:na+nc+ng;
vd=na+nc+ng+1:na+nc+ng+nd;
vb=na+nc+ng+nd+1:na+nc+ng+nd+nb;

PwNl=repmat(Pw,Nl,1);
PwNl=PwNl(ka1==1);

     v1=zeros(na,1);
     v1([RC.ACr;RC.AGr])=1;
     r=find(v1==1);
 
% aw1=sum(SCw(vc).*dV(vc));

 if isempty(RC.Cr)==0 || isempty(RC.Gr)==0
     [Bwc,Bwg,Blc,Blg,Sw(vc),Sw(vg),ndt,Q1,Q2,Qm,~]=fun1(RC,Pi,Sw,Cp,PR,TC,TG,A2C,...
         A2G,WonC,WonG,Uf,CpW,Pw,dt,dV,CR_cr,Qz,Qf,ndt0,Pi0,L,Lc,Lg,Ke);
     %  [Bc,Bg,SCw(vc),SCw(vg),ndt,Q1,Q2,Qm,dSS]=fun2(RC,Pi,SCw,SCp,PR,TC,TG,A2C,...
     %     A2G,WonC,WonG,Uf,CpW,Pw,dt,dV,CR_cr,Qz,Qf,ndt0,Pi0,L,Lc,Lg,Ke,CR,TM);
     
     
     b1wm=sparse(Won,ones(1,size(Won,1)),-W1.*PwNl,na,1);
     b1wd=sparse(WoD(:,1),ones(1,size(WoD,1)),-W1D.*PwNl,na,1);
     
     W2M=sparse(WonM,Won,W1,nw,na);
     W2D=sparse(WoD(:,3),WoD(:,1),W1D,nw,nd);
     W2B=sparse(nw,nb);

     b1wm=b1wm.*(sum(W2M(Qf==0,:),1)~=0)';
     b1wd=b1wd.*(sum(W2M(Qf==0,:),1)~=0)';
     
     WM1=[W2M,W2D,W2B];
     WM2=WM1';
     W3vec=sparse(WonM,1,W1,nw,1)+sparse(WoD(:,3),1,W1D,nw,1);
     WM3=-sparse(1:nw,1:nw,W3vec,nw,nw);
  
     Qf=Qzm1<0;
            
     bl=[b1wm,b1wd,b1wb]'+BLGY_GIM;
     Bl=bl'-sparse(r,ones(sum(v1),1),Blc,na+nd+nb,1)/dt-sparse(r,ones(sum(v1),1),Blg,na+nd+nb,1)/dt;
     
     A1=TL-sparse(Won,Won,W1,na,na)-sparse(1:na,1:na,Clp(va)+sum(b1gm(:,1:2))+sum(A2DL,2),na,na);
     D1=DL-sparse(1:nd,1:nd,sum(A2DL,1)+Clp(vd)'+sum(b1gd(:,1:2)),nd,nd)-sparse(WoD(:,1),WoD(:,1),W1D,nd,nd);                       %Матрица коэф. для двойной пор.
     B1=BB-sparse(1:nb,1:nb,sum(A2BL,1)+Clp(vb)'+b1gb,nb,nb);
        
     AM=[A1,   A2DL, A2BL;
        A2DL',  D1,  D2BL;
        A2BL', D2BL', B1];
     
     W2M=W2M(Qf~=0,:);
     WM2=WM2(:,Qf~=0);
     WM3=WM3(Qf~=0,Qf~=0);
     Qzm1(CR_cr.wn)=(Qm(:,1)+Qm(:,2))/dt;
     Pt=[Bl',Qzm1(Qf~=0)']/[AM,WM2;W2M,WM3];
     
     Pi(va)=Pt(va)';
     Pi(vd)=Pt(vd)';
     Pi(vb)=Pt(vb)';
     PwNl(Qf~=0)=Pt(na+nd+nb+1:end);
     
     %временно
     Bpc=zeros(size(Bwc));
     Bpg=zeros(size(Bwc));
    
 else
     Bwc=zeros(nc,1);
     Bwg=zeros(ng,1);
     
     Bpc=zeros(size(Bwc));
     Bpg=zeros(size(Bwc));
    
     ndt=1;
     Q1=zeros(size(WonC(:,3),1),5);
     Q2=zeros(size(WonG(:,3),1),5);
     Qm=zeros(0,5);
 end;

b_A2B=Soed2B(A2BW,A2BP,Pi,na,nb,va,vb);     % Связь пор с граничной областью
b_D2B=Soed2B(D2BW,D2BP,Pi,nd,nb,vd,vb);     % Связь трещин с граничной областью

     bwm=sparse(Won,ones(1,size(Won,1)),-W6.*(Pi(Won)-PwNl),na,1)...
         -b1gm(:,3).*(Pi(va)-GYData.GY_Pxy)-b1gm(:,4).*(Pi(va)-GYData.GY_Pz)-b_A2B(:,1);
     
     bwd=sparse(WoD(:,1),ones(1,size(WoD,1)),-W6D.*(Pi(WoD(:,1))-PwNl(WoD(:,3))),nd,1)...
         -b1gd(:,3).*(Pi(vd)-GYData.GY_Pxy(vad))-b1gd(:,4).*(Pi(vd)-GYData.GY_Pz(vad))-b_D2B(:,1);
    
     bpm=sparse(Won,ones(1,size(Won,1)),-W7.*(Pi(Won)-PwNl),na,1)...
         -b1gm(:,5).*(Pi(va)-GYData.GY_Pxy)-b1gm(:,6).*(Pi(va)-GYData.GY_Pz)-b_A2B(:,2);
     
     bpd=sparse(WoD(:,1),ones(1,size(WoD,1)),-W7D.*(Pi(WoD(:,1))-PwNl(WoD(:,3))),nd,1)...
         -b1gd(:,5).*(Pi(vd)-GYData.GY_Pxy(vad))-b1gd(:,6).*(Pi(vd)-GYData.GY_Pz(vad))-b_D2B(:,2);
     
     bw=[bwm;bwd];
     bp=[bpm;bpd];

     Bw=bw+sparse(r,ones(sum(v1),1),Bwc,na+nd,1)/dt+sparse(r,ones(sum(v1),1),Bwg,na+nd,1)/dt-Cwp([va,vd]).*(Pi([va,vd])-Pi0([va,vd]));
     Bp=bp+sparse(r,ones(sum(v1),1),Bpc,na+nd,1)/dt+sparse(r,ones(sum(v1),1),Bpg,na+nd,1)/dt-Cp([va,vd]).*Cwp([va,vd]).*(Pi([va,vd])-Pi0([va,vd]));
     tmp=sum(Bw);    

     AM2=[TW-sparse(1:na,1:na,sum(A2DW,2)),A2DW; A2DW',DW-sparse(1:nd,1:nd,sum(A2DW,1))];
     AM3=[TP,A2DP;A2DP',DP];
     
     Sw_old=Sw([va,vd]);
     Sw([va,vd])=Sw([va,vd])+dt*(AM2*Pi([va,vd])+Bw)./Cws([va,vd]);

     Cp([va,vd])=Sw_old.*Cp([va,vd])+dt*(AM3*Pi([va,vd])+Bp)./Cws([va,vd]);
     
     vad=[va,vd];
     v0=Sw(vad)==0;
     Cp(vad(v0==0))=Cp(vad(v0==0))./Sw(vad(v0==0));

    % Sw=Sw.*(Sw>=PR.aw(4)).*(Sw<=1-PR.aw(5))+(1-PR.aw(5)).*(Sw>1-PR.aw(5))+PR.aw(4).*(Sw<PR.aw(4));
    Sw=Sw.*(Sw>=0).*(Sw<=1)+(Sw>1); 
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