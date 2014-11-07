%%Ãîòîâèì ÐÈÃÈÑ
load('RIGIS','NB')
load('WellCoord','WXY')

NB(:,1)=N2C(NB(:,1));
NB(:,1)=PrepText(NB(:,1));
NB(:,1)=PoiskTire(NB(:,1));

NB(:,3)=N2C(NB(:,3));
NB(:,3)=PrepText(NB(:,3));
NB(:,3)=PoiskTire(NB(:,3));

LName=unique(NB(:,3));

for i=10:12
    im=strcmp(LName(i),NB(:,3));
    NB(im==1,:)=[];
end;


K=zeros(size(wn,1),9);
Sw=zeros(size(wn,1),9);
Mp=zeros(size(wn,1),9);
Hk=zeros(size(wn,1),9);
H=zeros(size(wn,1),9);

for i=1:size(wn,1)
  r=strcmp(wn(i),NB(:,1));
  nb=NB(r==1,:);
  
  for j=1:size(nb,1)
    for k=1:10  
     date=nb{j,k+3};
     img=isempty(date);
     date(img==1)=NaN;
     nb(j,k+3)={date};
    end;
  end;
  
  
    for l=1:9
        im=strcmp(LName(l),nb(:,3));
        if sum(im)>0
            K(i,l)=nb{im==1,12};
            Sw(i,l)=1-nb{im==1,9};
            Mp(i,l)=nb{im==1,10};
            Hk(i,l)=nb{im==1,7}-nb{im==1,6};
            Z(i,l)=(nb{im==1,5}+nb{im==1,4})/2;
            H(i,l)=nb{im==1,5}-nb{im==1,4};
        end;
    end;
end;

[K,Sw,Mp,Hk,H]=filter1(K,Sw,Mp,Hk,H,WXY);