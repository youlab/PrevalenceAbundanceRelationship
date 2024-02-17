clear;
clc;
close all;
load('SevenDiseases.mat');
load("AMRIndex.mat");
CC=linspecer(3);
AMR_fitness=zeros(2,length(GeneID));
AMR_prev=zeros(2,length(GeneID));
AMR_abund=zeros(2,length(GeneID));

Data=CrohnsDisease_None;
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
p=polyfit(A,B,1)
subplot(2,1,1);
plot(A,B,'r.','markersize',3,'color',CC(1,:));hold on;
set(gca,'fontsize',10);

Data=CrohnsDisease;
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
p=polyfit(A,B,1)
subplot(2,1,2);
plot(A,B,'g.','markersize',3,'color',CC(2,:));hold on;
set(gca,'fontsize',10);

set(gcf,'position',[100 100 200 450]);

saveas(gcf,'CrohnsDisease_PowerLaw.fig');
saveas(gcf,'CrohnsDisease_PowerLaw.pdf');
