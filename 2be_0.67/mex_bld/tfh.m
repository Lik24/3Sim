function C=tfh(A,B)
A1=sparse(A);
C=A1*B;
C=full(C);