clear;
close all;
clc;
load('D:\GeneCorrelationAnalysis\HumanSkin\RawData\RawDataGeneNose.mat');
Data0=SampleHumanSkin;
Data0=Data0./sum(Data0,1);

NNs=10:5:size(Data0,2);
repeat=20;
Slope=zeros(length(NNs),repeat);

for i=1:length(NNs)
    NN=NNs(i);
    NN
    for dfg=1:repeat
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
        p=polyfit(x,y,1);
        Slope(i,dfg)=p(1);
    end
end
C=linspecer(2);
errorbar(NNs,mean(Slope,2),std(Slope,0,2),'o','markersize',5,'capsize',0,'color',C(1,:));hold on;

set(gca,'fontsize',8);
set(gcf,'position',[100 100 220 220]);
xlim([5,size(Data0,2)]);
xlabel('sample number','fontsize',20);
ylabel('\theta value','fontsize',20);
saveas(gcf,'HumanNose.fig');
saveas(gcf,'HumanNose.pdf');
