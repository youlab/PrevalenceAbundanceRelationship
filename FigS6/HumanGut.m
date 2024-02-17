clear;
close all;
clc;
load('D:\GeneCorrelationAnalysis\HumanFourCountries\RawData\RawDataGene.mat');

Data0=SampleAll;
Data0=Data0./sum(Data0,1);
Prev0=0;
Abun0=0;
for j=1:size(Data0,1)
    Prev0(j)=nnz(Data0(j,:))./size(Data0,2);
    Abun0(j)=mean(nonzeros(Data0(j,:)));
    Phi0(j)=-log(1-Prev0(j))./Prev0(j);
end
x=log(Phi0(Prev0>0&Prev0<1));
y=log(Abun0(Prev0>0&Prev0<1));
mm0=mean(y-x);

repeat=100;
SubSample=20;
Prev=zeros(size(Data0,1),repeat);
Abun=zeros(size(Data0,1),repeat);
Phi=zeros(size(Data0,1),repeat);
for i=1:repeat
    temp=randperm(size(Data0,2));
    index=temp(1:SubSample);
    Data=Data0(:,index);
    for j=1:size(Data,1)
        Prev(j,i)=nnz(Data(j,:))./size(Data,2);
        Abun(j,i)=mean(nonzeros(Data(j,:)));
        Phi(j,i)=-log(1-Prev(j,i))./Prev(j,i);
    end
    PPrev=Prev(:,i);
    AAbun=Abun(:,i);
    PPhi=Phi(:,i);
    x=log(PPhi(PPrev>0&PPrev<1));
    y=log(AAbun(PPrev>0&PPrev<1));
    mm(i)=mean(y-x);
end

GeneName=Name;
GeneAbundance=SampleAll;
clear Name;
clear SampleAll;
TxtFile=fileread('D:\GeneCorrelationAnalysis\HumanFourCountries\RawData\KEGGMapper_Brite.txt');


BriteTarget='<div id="object';
BriteStrIndex=strfind(TxtFile,BriteTarget);

for i=1:length(BriteStrIndex)
    Name(i,:)=TxtFile([BriteStrIndex(i)+length(BriteTarget):BriteStrIndex(i)+length(BriteTarget)+6]); 
end
Name=table(Name);
BriteStrIndex(length(BriteStrIndex)+1)=length(TxtFile);
SampleAll=0*ones(size(Name,1),size(GeneAbundance,2));


Fitness0=zeros(length(BriteStrIndex)-1,size(GeneName,1));
Fitness=zeros(length(BriteStrIndex)-1,size(GeneName,1),repeat);

 for i=1:size(GeneName,1)
     
    GeneStrIndex=strfind(TxtFile,[char(GeneName{i,:}) '</a></dt>']);
    for j=1:length(GeneStrIndex)
        for k=1:length(BriteStrIndex)-1
            if GeneStrIndex(j)>BriteStrIndex(k)&&GeneStrIndex(j)<BriteStrIndex(k+1)
                Fitness0(k,i)=log(Abun0(i))-log(Phi0(i))-mm0;
                for fgh=1:repeat
                    Fitness(k,i,fgh)=log(Abun(i,fgh))-log(Phi(i,fgh))-mm(fgh);
                end
                break;
            end
        end
    end
 end

for i=1:length(Fitness0(:))
    if Fitness0(i)==-Inf||Fitness0(i)==0
            Fitness0(i)=NaN;    
    end
end

for i=1:length(Fitness(:))
    if Fitness(i)==-Inf||Fitness(i)==0
            Fitness(i)=NaN;    
    end
end

C=linspecer(3);
 ii=42;
 figure(1);
 for i=1:size(Fitness0,2)
    if ~isnan(Fitness0(ii,i))
        i
        errorbar(Fitness0(ii,i),mean(Fitness(ii,i,:),'omitnan'),std(Fitness(ii,i,:),'omitnan'),'o','markersize',10,'CapSize',0);hold on;
    end
 end

 set(gca,'fontsize',12);
 xlabel('\chi value in the full dataset','fontsize',14);
 ylabel('\chi values in the subsampling','fontsize',14);
 set(gcf,'position',[100 100 250 250]);
saveas(gcf,'HumanGutAMR.fig');
saveas(gcf,'HumanGutAMR.pdf');


 ii=14;
figure(2);
 for i=1:size(Fitness0,2)
    if ~isnan(Fitness0(ii,i))
        i
        errorbar(Fitness0(ii,i),mean(Fitness(ii,i,:),'omitnan'),std(Fitness(ii,i,:),'omitnan'),'o','markersize',10,'CapSize',0);hold on;
    end
 end

 set(gca,'fontsize',12);
 xlabel('\chi value in the full dataset','fontsize',14);
 ylabel('\chi values in the subsampling','fontsize',14);
 set(gcf,'position',[100 100 250 250]);
saveas(gcf,'HumanGutPhoto.fig');
saveas(gcf,'HumanGutPhoto.pdf');