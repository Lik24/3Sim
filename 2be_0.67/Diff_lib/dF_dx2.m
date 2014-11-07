function [dF ] = dF_dx2(dTp,dW,Pw,P,Won,na,D,flag1 )

if flag1==1 %dF1_dx2
    dF=dTp+sparse(Won,Won,dW.*(Pw-P),na,na);
else        %dF2_dx2
    dF=D+dTp+sparse(Won,Won,dW.*(Pw-P),na,na);
end

end

