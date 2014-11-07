function [td,tz,tim_step,NDay_D,NDay_Z]=History_Time(DOB,ZAK)
Day=[31,28,31,30,31,30,31,31,30,31,30,31];

day_d=zeros(size(DOB,1),1);
mes_d=zeros(size(DOB,1),1);
god_d=zeros(size(DOB,1),1);

day_z=zeros(size(ZAK,1),1);
mes_z=zeros(size(ZAK,1),1);
god_z=zeros(size(ZAK,1),1);

for i=1:size(DOB,1)
tex=DOB{i};    
day_d(i,1)=str2num(tex(1:2));
mes_d(i,1)=str2num(tex(4:5));
god_d(i,1)=str2num(tex(7:10));
end

for i=1:size(ZAK,1)
tex=ZAK{i};    
day_z(i,1)=str2num(tex(1:2));
mes_z(i,1)=str2num(tex(4:5));
god_z(i,1)=str2num(tex(7:10));
end

god0=min([god_d;god_z]);
god_d=god_d-god0;
god_z=god_z-god0;

for i=1:size(day_d,1)
    td(i)=day_d(i)+sum(Day(1:mes_d(i)-1))+god_d(i)*365;
    NDay_D(i,1)=Day(mes_d(i));
    NDay_D(i,2)=day_d(i);
end;

for i=1:size(day_z,1)
    tz(i)=day_z(i)+sum(Day(1:mes_z(i)-1))+god_z(i)*365;
    NDay_Z(i,1)=Day(mes_z(i));
    NDay_Z(i,2)=day_z(i);
end;
min(td)
min(tz)
mtzd=min([td,tz]);
td=td-mtzd+1;
tz=tz-mtzd+1;

tim_step=max([td,tz])-min([td,tz]);