%%Готовим историю
load('DOB_ZAK_in','WN','DOB','ZAK')
load('DOB_ZAK_in','DOBM','ZAKM')
wn=N2C(WN);

wn=N2C(wn);
wn_dob=N2C(DOB(:,1));
wn_zak=N2C(ZAK(:,1));

wn=PrepText(wn);
wn_dob=PrepText(wn_dob);
wn_zak=PrepText(wn_zak);
DOB(:,2)=PrepText(DOB(:,2));
ZAK(:,2)=PrepText(ZAK(:,2));

% for i=1:size(DOB,1)
%  for j=1:9   
%   if isempty(DOB{i,j+3})==0  
%    DOBM(i,j)=DOB{i,j+3};
%   else
%    DOBM(i,j)=nan;
%   end
%  end;
% end;
% 
% for i=1:size(ZAK,1)
%  for j=1:7   
%   if isempty(ZAK{i,j+3})==0  
%    ZAKM(i,j)=ZAK{i,j+3};
%   else
%    ZAKM(i,j)=nan;
%   end
%  end;
% end;

LName=unique([DOB(:,2);ZAK(:,2)]);

r1=strcmp('K',DOB(:,2));
r2=strcmp('K',ZAK(:,2));

[DB,ZK,tim_step,NDay_D,NDay_Z]=History_Time(DOB(r1,3),ZAK(r2,3));

dob=DOB(r1,:);
dobM=DOBM(r1,:);
wn_dob=wn_dob(r1);

zak=ZAK(r2,:);
zakM=ZAKM(r2,:);
wn_zak=wn_zak(r2);

dob_oil=zeros(size(wn,1),tim_step);
dob_lik=zeros(size(wn,1),tim_step);
dob_gas=zeros(size(wn,1),tim_step);
dob_Ppl=zeros(size(wn,1),tim_step);
dob_Pw=zeros(size(wn,1),tim_step);
dob_time=zeros(size(wn,1),tim_step);
dob_time2=zeros(size(wn,1),tim_step);

dob_oil2=zeros(size(wn,1),tim_step);
dob_lik2=zeros(size(wn,1),tim_step);
dob_gas2=zeros(size(wn,1),tim_step);
dob_Ppl2=zeros(size(wn,1),tim_step);
dob_Pw2=zeros(size(wn,1),tim_step);

zak_lik=zeros(size(wn,1),tim_step);
zak_Ppl=zeros(size(wn,1),tim_step);
zak_Pw=zeros(size(wn,1),tim_step);
zak_Pbuf=zeros(size(wn,1),tim_step);
zak_time=zeros(size(wn,1),tim_step);
zak_time2=zeros(size(wn,1),tim_step);

zak_lik2=zeros(size(wn,1),tim_step);
zak_Ppl2=zeros(size(wn,1),tim_step);
zak_Pw2=zeros(size(wn,1),tim_step);
zak_Pbuf2=zeros(size(wn,1),tim_step);

for i=1:size(wn,1)
  r1=strcmp(wn(i),wn_dob);
  dob=dobM(r1,:);
  
  r2=strcmp(wn(i),wn_zak);
  zak=zakM(r2,:);
  
  if isempty(dob(:,3))==0
      dtime=DB(r1);
      time=dtime;
      
      Dob_data=dob;
      dob_oil(i,time)=Dob_data(:,1);
      dob_lik(i,time)=Dob_data(:,3);
      dob_gas(i,time)=Dob_data(:,4);
      dob_Ppl(i,time)=Dob_data(:,7);
      dob_Pw(i,time)=Dob_data(:,8);
      dob_time(i,time)=NDay_D(r1,1);
      dob_time2(i,time)=NDay_D(r1,2);
  end;
  
  if isempty(zak(:,3))==0
      dtime=ZK(r2);
      time=dtime;
      
      Zak_data=zak;
      zak_lik(i,time)=Zak_data(:,1);
      zak_Ppl(i,time)=Zak_data(:,4);
      zak_Pw(i,time)=Zak_data(:,5);
      zak_Pbuf(i,time)=Zak_data(:,6);
      zak_time(i,time)=NDay_Z(r2,1);
      zak_time2(i,time)=NDay_Z(r2,2);
  end;
end;

for i=1:size(wn,1)
  A=dob_time(i,:);
  A1=dob_time2(i,:);  
  r=find(A1);
  for j=1:size(r,2)
    dob_oil2(i,[r(j)-A1(r(j))+1:r(j)+A(r(j))])=dob_oil(i,r(j))/A(r(j));
    dob_lik2(i,[r(j)-A1(r(j))+1:r(j)+A(r(j))])=dob_lik(i,r(j))/A(r(j));
    dob_gas2(i,[r(j)-A1(r(j))+1:r(j)+A(r(j))])=dob_gas(i,r(j))/A(r(j));
    dob_Ppl2(i,[r(j)-A1(r(j))+1:r(j)+A(r(j))])=dob_Ppl(i,r(j));
    dob_Pw2(i,[r(j)-A1(r(j))+1:r(j)+A(r(j))])=dob_Pw(i,r(j));
  end;
  
  A=zak_time(i,:);
  A1=zak_time2(i,:);  
  r=find(A1);
  for j=1:size(r,2)
    zak_lik2(i,[r(j)-A1(r(j))+1:r(j)+A(r(j))])=zak_lik(i,r(j))/A(r(j));
    zak_Pbuf2(i,[r(j)-A1(r(j))+1:r(j)+A(r(j))])=zak_Pbuf(i,r(j));
    zak_Ppl2(i,[r(j)-A1(r(j))+1:r(j)+A(r(j))])=zak_Ppl(i,r(j));
    zak_Pw2(i,[r(j)-A1(r(j))+1:r(j)+A(r(j))])=zak_Pw(i,r(j));
  end;
end;

indate=sum([dob_oil2; dob_lik2; dob_gas2; dob_Ppl2; dob_Pw2; zak_lik2; zak_Pbuf2; zak_Ppl2; zak_Pw2]);
ri=find(indate~=0);

dob_oil2=dob_oil2(:,ri(1):end);
dob_lik2=dob_lik2(:,ri(1):end);
dob_gas2=dob_gas2(:,ri(1):end);
dob_Ppl2=dob_Ppl2(:,ri(1):end);
dob_Pw2=dob_Pw2(:,ri(1):end);

zak_lik2=zak_lik2(:,ri(1):end);
zak_Pbuf2=zak_Pbuf2(:,ri(1):end);
zak_Ppl2=zak_Ppl2(:,ri(1):end);
zak_Pw2=zak_Pw2(:,ri(1):end);

