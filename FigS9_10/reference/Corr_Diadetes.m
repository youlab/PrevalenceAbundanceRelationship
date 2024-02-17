clear;
clc;
close all;
load('SevenDiseases.mat');
load("AMRIndex.mat");
AMR_fitness=zeros(2,length(GeneID));
AMR_prev=zeros(2,length(GeneID));
AMR_abund=zeros(2,length(GeneID));

Data=Type2Diadetes_None;
Data=Data./sum(Data,1);
prev=0;
abund=0;
phi=0;
for k=1:size(Data,1)
    prev(k)=nnz(Data(k,:))/size(Data,2);
    abund(k)=mean((Data(k,:)));
end
phi=-log(1-prev);
A=log(phi((prev<1)&(prev>0)));
B=log(abund((prev<1)&(prev>0)));
fitness=log(abund)-log(phi)-mean(B-A);
AMR_fitness(1,:)=fitness;
AMR_prev(1,:)=prev;
AMR_abund(1,:)=abund;

Data=Type2Diadetes;
Data=Data./sum(Data,1);
prev=0;
abund=0;
phi=0;
for k=1:size(Data,1)
    prev(k)=nnz(Data(k,:))/size(Data,2);
    abund(k)=mean((Data(k,:)));
end
phi=-log(1-prev);
A=log(phi((prev<1)&(prev>0)));
B=log(abund((prev<1)&(prev>0)));
fitness=log(abund)-log(phi)-mean(B-A);
AMR_fitness(2,:)=fitness;
AMR_prev(2,:)=prev;
AMR_abund(2,:)=abund;

plot(AMR_fitness(1,:),AMR_fitness(2,:),'k.');hold on;
yy=-6:0.1:4;
plot(0*yy,yy,'r--');hold on;
xx=-5:0.1:5;
plot(xx,0*xx,'r--');hold on;

set(gcf,'position',[100 100 200 200]);
saveas(gcf,'Corr_Diabetes.fig');
saveas(gcf,'Corr_Diabetes.pdf');