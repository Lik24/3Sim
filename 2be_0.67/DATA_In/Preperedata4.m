%%Ãîòîâèì ÐÈÃÈÑ
load('RIGIS','NB','NBS1')
load('WCoord','WXY')

NBS1(:,1)=N2C(NBS1(:,1));
NBS1(:,1)=PrepText(NBS1(:,1));
NBS1(:,1)=PoiskTire(NBS1(:,1));

NBS1(:,2)=N2C(NBS1(:,2));
NBS1(:,2)=PrepText(NBS1(:,2));
NBS1(:,2)=PoiskTire(NBS1(:,2));

NB(:,3)=N2C(NB(:,3));
NB(:,3)=PrepText(NB(:,3));
NB(:,3)=PoiskTire(NB(:,3));

LName=unique(NB(:,3));

for i=10:12
    im=strcmp(LName(i),NBS1(:,2));
    NBS1(im==1,:)=[];
end;


% K=zeros(size(wn,1),9);
% Sw=zeros(size(wn,1),9);
% Mp=zeros(size(wn,1),9);
% Hk=zeros(size(wn,1),9);
% H=zeros(size(wn,1),9);

for i=1:size(wn,1)
  r=strcmp(wn(i),NBS1(:,1));
  nb=NBS1(r==1,:);
  
%   for j=1:size(nb,1)
%     for k=1:10  
%      date=nb{j,k+3};
%      img=isempty(date);
%      date(img==1)=NaN;
%      nb(j,k+3)={date};
%     end;
%   end;
  
  
    for l=1:9
        im=strcmp(LName(l),nb(:,2));
        if sum(im)>0
            nb{im==1,10}
            perf(i,l)=nb{im==1,10};
            Sw(i,l)=1-nb{im==1,11};
            Mp(i,l)=nb{im==1,10};
            Hk(i,l)=nb{im==1,7}-nb{im==1,6};
            H(i,l)=nb{im==1,5}-nb{im==1,4};
        end;
    end;
end;

[K,Sw,Mp,Hk,H]=filter1(K,Sw,Mp,Hk,H,WXY);