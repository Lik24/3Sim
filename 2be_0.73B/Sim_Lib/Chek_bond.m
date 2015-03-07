function [fl_out,Pw,Q,qf]=Chek_bond(Pw,Ppl,W1,uf,qf,PQ)
Q=W1.*(Ppl'-Pw);

fl(:,1)=Pw<PQ(1,1);
fl(:,2)=Pw>PQ(1,2);
fl(:,3)=(Q<PQ(2,1)).*(uf==1);
fl(:,4)=(Q>PQ(2,2)).*(uf==1);
fl(:,5)=(Q<PQ(3,1)).*(uf==-1);
fl(:,6)=(Q>PQ(3,2)).*(uf==-1);

for i=1:size(uf,1)
    if uf(i)==1
        if qf(i)==0
            Q(i)=PQ(2,1)*(fl(i,3)==1)+Q(i)*(fl(i,3)==0);
            Q(i)=PQ(2,2)*(fl(i,4)==1)+Q(i)*(fl(i,4)==0);
            qf(i)=(fl(i,3)+fl(i,4))~=0;
        else
            Pw(i)=PQ(1,1)*(fl(i,1)==1)+Pw(i)*(fl(i,1)==0);
            Pw(i)=PQ(1,2)*(fl(i,2)==1)+Pw(i)*(fl(i,2)==0);
            qf(i)=(fl(i,1)+fl(i,2))==0;
        end
    else
        if qf(i)==0
            Q(i)=PQ(3,1)*(fl(i,5)==1)+Q(i)*(fl(i,5)==0);
            Q(i)=PQ(3,2)*(fl(i,6)==1)+Q(i)*(fl(i,6)==0);
            qf(i)=(fl(i,5)+fl(i,6))~=0;
        else
            Pw(i)=PQ(1,1)*(fl(i,1)==1)+Pw(i)*(fl(i,1)==0);
            Pw(i)=PQ(1,2)*(fl(i,2)==1)+Pw(i)*(fl(i,2)==0);
            qf(i)=(fl(i,1)+fl(i,2))==0;
        end
    end;
end;

fl_out=sum(fl(:))>0;
