function [TL,TW,TG]=Potok_MKT_GPU(Lm,Bm,Ke,P,Sw,H,r,c,as,aw)
n=size(Lm,1);
L=full(Lm(r+(c-1)*n));
B=full(Bm(r+(c-1)*n));
H=full(Hm(r+(c-1)*n));

gL=gpuArray(L);
gB=gpuArray(B);
gSw=gpuArray(Sw);
gP=gpuArray(P);
gK=gpuArray(Ke);
gH=gpuArray(H);

vP=gP(r)>=gP(c);
Swc=gSw(r);
Swl=gSw(c);

gSwe=Swc.*vP+Swl.*(vP==0);
Kfo=Sat_cal(gSwe,2,1,as,aw); %oil
Kfw=Sat_cal(gSwe,1,1,as,aw); %water

T=gK.*gB.*gH./gL*2;

gTw=T.*Kfw;
gTo=T.*Kfo;

Tw=gather(gTw);
To=gather(gTo);

Tw1=sparse((r+(c-1)*n),ones(size(r)),Tw,n*n,1);
To1=sparse((r+(c-1)*n),ones(size(r)),To,n*n,1);

Tw=reshape(Tw1,n,n);
To=reshape(To1,n,n);

TL=To+Tw-sparse(1:n,1:n,sum(To+Tw),n,n);
TW=Tw-sparse(1:n,1:n,sum(Tw),n,n);
TG=1;