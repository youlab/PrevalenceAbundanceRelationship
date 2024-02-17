clear;
clc;
load('AMRIndex.mat');
load('SevenDiseases.mat');
AMR_fitness=zeros(2,length(AMRIndex));
AMR_prev=zeros(2,length(AMRIndex));
AMR_abund=zeros(2,length(AMRIndex));

for i=1:2
    if i==1
        Data=CrohnsDisease_None;
    end
    if i==2
        Data=CrohnsDisease;
    end
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
    AMR_fitness(i,:)=fitness(AMRIndex);
    AMR_prev(i,:)=prev(AMRIndex);
    AMR_abund(i,:)=abund(AMRIndex);
end

for i=1:length(AMR_fitness(:))
    if AMR_fitness(i)==-Inf
        AMR_fitness(i)=NaN;
    end
end
subplot(3,1,1);
boxplot(AMR_prev','symbol','');
set(gca,'fontsize',10);
ylim([0.7,1.05]);
xticklabels([]);

subplot(3,1,2);
boxplot(AMR_abund','symbol','');
set(gca,'fontsize',10);
ylim([-1,4]*10^(-4));
xticklabels([]);

subplot(3,1,3);
boxplot(AMR_fitness','symbol','');hold on;
plot([-1:2],0*[-1:3],'k--');
set(gca,'fontsize',10);

set(gcf,'position',[100 100 200 500]);

saveas(gcf,'AMR_CrohnsDisease.fig');
saveas(gcf,'AMR_CrohnsDisease.pdf');


