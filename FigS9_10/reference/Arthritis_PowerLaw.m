clear;
clc;
close all;
load('SevenDiseases.mat');
load("AMRIndex.mat");
CC=linspecer(3);
AMR_fitness=zeros(3,length(GeneID));
AMR_prev=zeros(3,length(GeneID));
AMR_abund=zeros(3,length(GeneID));

Data=Arthritis_None;
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
subplot(3,1,1);
plot(A,B,'r.','markersize',3,'color',CC(1,:));hold on;
set(gca,'fontsize',10);

Data=Arthritis_Moderate;
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
subplot(3,1,2);
plot(A,B,'g.','markersize',3,'color',CC(2,:));hold on;
set(gca,'fontsize',10);

Data=Arthritis_High;
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
subplot(3,1,3);
plot(A,B,'b.','markersize',3,'color',CC(3,:));hold on;
set(gca,'fontsize',10);

set(gcf,'position',[100 100 200 700]);

saveas(gcf,'Arthritis_PowerLaw.fig');
saveas(gcf,'Arthritis_PowerLaw.pdf');
