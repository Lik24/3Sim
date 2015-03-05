function [FN,Fc,dS,Sc]=simbol2mat(f,F,FF)

syms Sw
ko=(1-Sw).^2;
kw=0.1*Sw.^2;
mu=[1,300];
mp=0.34;
K2=0.65;
Swr=0.25;
Sor=1-Swr-K2*(1-Swr);
gam=mu(1)/mu(2);
H=126.8895*(1-0.8012);
V0=mp*H*A*(1-Swr);

f=kw/(gam*ko+kw);
F=diff(f,Sw);
FF=diff(F,Sw);

dS=0:0.005:1;
fn=(subs(f,dS)-subs(f,0))./(dS-0);
fm=max(eval(fn));%
[r,c]=find(fm==fn);
Sc=dS(c);
Fc=eval(subs(F,Sc));

FN(:,1)=eval(subs(f,dS));
FN(:,2)=eval(subs(F,dS));
FN(:,3)=eval(subs(FF,dS));
