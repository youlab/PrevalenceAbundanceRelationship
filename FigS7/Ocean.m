clear;
close all;
clc;
load('D:\GeneCorrelationAnalysis\TaraOcean\RawData\RawDataGene.mat');
Data0=TotalReads(2:end,:);
Data0=Data0./sum(Data0,1);

NNs=10:10:size(Data0,2);
repeat=20;
Chi=zeros(length(NNs),size(Data0,1));

for i=1:length(NNs)
    NN=NNs(i);
    NN
    temp=randperm(size(Data0,2));
    index=temp(1:NN);
    Data=Data0(:,index);
    Prev=0;
    Abun=0;
    for j=1:size(Data,1)
        Prev(j)=nnz(Data(j,:))./size(Data,2);
        Abun(j)=mean((Data(j,:)));
    end
    Phi=-log(1-Prev);
    x=log(Phi(Prev>0&Prev<1));
    y=log(Abun(Prev>0&Prev<1));
    Chi(i,:)=log(Abun)-log(Phi)-mean(y-x);
end
C=linspecer(3);

for i=1:length(Chi(:))
    if Chi(i)==-Inf||Chi(i)==0
        Chi(i)=NaN;
    end
end
plot(NNs,Chi);hold on;
errorbar(NNs,mean(Chi,2,'omitnan'),std(Chi,0,2,'omitnan'),'ko-','markersize',10,'linewidth',2,'capsize',0);hold on;
set(gca,'fontsize',10);
xlabel('Number of samples','FontSize',12);
ylabel('\chi values','fontsize',12);
xlim([5 size(Data0,2)]);
ylim([-6,6]);
set(gcf,'position',[100 100 250 250]);
saveas(gcf,'Ocean.fig');
saveas(gcf,'Ocean.pdf');


