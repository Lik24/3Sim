load('history_deb_day','dob_lik2','zak_lik2');

Qlz(:,:)=Q(:,2,:)+Q(:,1,:);
sQlz=cumsum(Qlz');

Qlzf=dob_lik2-zak_lik2;
sQlzf=cumsum(Qlzf');
T=1:1035;  

for i=33:-1:1

uf=WData.Uf(i,:);
 
figure(i);

%subplot(2,1,1),plot(T(uf~=0),Pw1(i,uf~=0),T(uf~=0),Pw(i,uf~=0),T(uf~=0),Pw_f(i,uf~=0))
subplot(2,1,1),plot(T(uf~=0),Pw3(i,uf~=0),T(uf~=0),Pw4(i,uf~=0),T(uf~=0),Pw_f(i,uf~=0))
subplot(2,1,2),plotyy([T',T'],[sQlz(:,i),sQlzf(:,i)],[T',T'],[Qlz(i,:)',Qlzf(i,:)'])
end

% for i=33:-1:1
% figure(i),plot(1:1035,Pw3(i,:),1:1035,Pw4(i,:),1:1035,Pw5(i,:),1:1035,Pw_f(i,:))
% end