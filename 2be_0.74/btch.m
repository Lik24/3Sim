fg=load('a_uni.mat');
a_uni=fg.a_uni50;
parfor i=1:1
 CD(:,i)=main1(a_uni(:,28));   
end