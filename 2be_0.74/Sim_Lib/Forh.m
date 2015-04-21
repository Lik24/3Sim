function [T] = Forh(T,kms, Ro, Kf, K, mu, dP)

F=forge(kms, Ro,Kf,K,mu, dP);
T=T.*(F.*(F~=0)+1.*(F==0)); 
end

function F=forge(kms, Ro, Kf, K, mul, dP)
 B=kms*Ro*(Kf.*K).^2./mul^2.*dP;
 F=ones(size(B));
 F(B~=0)=(((4.*B(B~=0)+1).^0.5-1)./(2.*B(B~=0)));
end

