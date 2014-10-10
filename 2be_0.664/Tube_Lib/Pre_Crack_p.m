function [CR_ind]=Pre_Crack_p(RC,na,T,Wf,Won,WonM,Wn,CR_GRUP,C,G,A2C,A2G,K,WonC,WonG)

na=RC.na;
nc=RC.nc;
vc=na+1:na+nc;

% ng=RC.ng;
[r4,c4]=find(A2G);

v1=zeros(na,1);

r=RC.Arc(:,1);
c=RC.Arc(:,2);

L_grup=unique(CR_GRUP(:,2));
soc=0;
for l=1:size(L_grup,1)
    l_grup=CR_GRUP(L_grup(l)==CR_GRUP(:,2),1);
    A2C_L=A2C(:,L_grup(l)==CR_GRUP(:,2));
    C_L=C(:,L_grup(l)==CR_GRUP(:,2));
    C_L=C_L(L_grup(l)==CR_GRUP(:,2),:);
    vc_L=vc(L_grup(l)==CR_GRUP(:,2));
    C_grup=unique(l_grup);
    for j=1:size(C_grup,1)
        soc=soc+1;
        A2C_CrS=A2C_L(:,C_grup(j)==l_grup);
        C_CrS=C_L(:,C_grup(j)==l_grup);
        C_CrS=C_CrS(C_grup(j)==l_grup,:);
        vc_CrS=vc_L(C_grup(j)==l_grup);
        
        A2C_cell(soc)={A2C_CrS};
        [r1,c1]=find(A2C_CrS);
        
        CR_rc=Narez_Index(v1,r1,na,r,c,T,Wf,Won,WonM,Wn,A2C_CrS,A2G);
        
        [r2,c2]=find(C_CrS>0);
        RC.Cr=r2;
        RC.Cc=c2;
        nc=size(C_CrS,1);
        RC2=sparse(r2,c2,1);
        U=triu(RC2);
        [r2_h,c2_h]=find(U);
        RC.Cr2=r2_h;
        RC.Cc2=c2_h;
               
        [r3,c3]=find(G>0);
        RC.Gr=r3;
        RC.Gc=c3;
        ng=size(G,1);
        RC3=sparse(r3,c3,1);
        U=triu(RC3);
        [r5_h,c5_h]=find(U);
        RC.Gr2=r5_h;
        RC.Gc2=c5_h;
        
        RC.ACr=r1;
        RC.ACc=c1;
        
        RC.AGr=r4;
        RC.AGc=c4;
        
        k=0;
        for i=1:size(WonC,1)
         ri=find(WonC(i,1)==vc_CrS);
         wonC(1,1:3)=1;
         wonC(1,:)=[];
         if isempty(ri)==0
           k=k+1;  
           wonC(k,:)=WonC(ri,:);
         end;
        end;
        
        CR_ind(soc,1)={CR_rc};
        CR_ind(soc,2)={C_CrS(r2_h+(c2_h-1)*nc)};
        CR_ind(soc,3)={G};
        CR_ind(soc,4)={A2C_CrS(r1+(c1-1)*na).*K(r1)};
        CR_ind(soc,5)={A2G};
        CR_ind(soc,6)={RC};
        CR_ind(soc,7)={vc_CrS};
        CR_ind(soc,8)={wonC};
        CR_ind(soc,9)={WonG};
        
    end;
end;
A2C_cell;


