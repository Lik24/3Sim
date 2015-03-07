function [D,A2D,dVd,pd,WonD,Ld]=DobPor(XY,d2,Nl,Won,r0,ka)
D=[];
A2D=ones(3*size(XY,1),1);
A2D(:,1)=[];
dVd=[];
pd=0;
WonD=0;
Ld=0;
