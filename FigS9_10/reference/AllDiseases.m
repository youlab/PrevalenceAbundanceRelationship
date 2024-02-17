clear;
clc;
load('AMRIndex.mat');
load('SevenDiseases.mat');
AMR_fitness=zeros(9,length(AMRIndex));
AMR_prev=zeros(9,length(AMRIndex));
AMR_abund=zeros(9,length(AMRIndex));

for i=1:9
    if i==10
        Data=Arthritis_None;
    end
    if i==11
        Data=ColorectalCancer_Controls;
    end
    if i==12
        Data=CrohnsDisease_None;
    end
    if i==13
        Data=LiverCirrhosis_None;
    end
    if i==1
        Data=Arthritis_Moderate;
    end
    if i==2
        Data=Arthritis_High;
    end
    if i==3
        Data=ColorectalCancer_Carcinoma;
    end
    if i==4
        Data=ColorectalCancer_AdvancedAdenoma;
    end
    if i==5
        Data=CrohnsDisease;
    end
    if i==6
        Data=LiverCirrhosis;
    end
    if i==7
        Data=Obesity;
    end
    if i==8
        Data=Type2Diadetes;
    end
    if i==9
        Data=UlcerativeColitis;
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
plot([-1:10],0*[-1:10],'k--');
set(gca,'fontsize',10);

set(gcf,'position',[100 100 400 500]);

saveas(gcf,'AllDiseases.fig');
saveas(gcf,'AllDiseases.pdf');


