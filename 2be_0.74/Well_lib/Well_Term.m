function [W1i,W1I,W1o,W2,W2T]=Well_Term(Q,TeW,Won,Ro,Cp,na,ndt,La,Wf,dt)

W1_in=Q(:,1).*Cp(1).*Ro(1).*TeW;
W1_In=Q(:,1).*Cp(1).*Ro(1);
W1_out=(Q(:,2)-Q(:,3)).*Cp(1).*Ro(1)+Q(:,3).*Cp(2).*Ro(2);
%W1_out=Q(:,2);

W1i=-sparse(Won,ones(size(Won,1),1),W1_in,na,1)/ndt;
W1I=-sparse(Won,ones(size(Won,1),1),W1_In,na,1)/ndt;
W1o=-sparse(Won,ones(size(Won,1),1),W1_out,na,1)/ndt;

W2_int=Wf*La(4)/0.5*dt/ndt;
W2=sparse(Won,ones(size(Won,1),1),W2_int,na,1);
W2T=sparse(Won,ones(size(Won,1),1),W2_int.*TeW,na,1);