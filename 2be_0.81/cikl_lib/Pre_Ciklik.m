function Psi=Pre_Ciklik(K,H,rc)
    
    h1=H(rc(:,2));
    h2=H(rc(:,1));
    
    k1=K(rc(:,2));
    k2=K(rc(:,1));
    
    Psi=(k1-k2).*h1.*h2./(k1.*h1+k2.*h2);
    Psi=sparse(rc(:,1),rc(:,2),Psi);