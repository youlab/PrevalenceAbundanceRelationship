clear;
clc;
load('RawDataGene.mat');
Data=SampleOcean;
Data=Data./sum(Data,1);
Prev=0;
Abun=0;
for j=1:size(Data,1)
    Prev(j)=nnz(Data(j,:))./size(Data,2);
    Abun(j)=mean((Data(j,:)));
    Phi(j)=-log(1-Prev(j));
end
x=log(Phi(Prev>0&Prev<1));
y=log(Abun(Prev>0&Prev<1));
mm=mean(y-x);

GeneName=Name;
GeneAbundance=Data;
clear Name;
clear SampleOcean;
TxtFile=fileread('KEGGMapper_Brite.txt');


BriteTarget='<div id="object';
BriteStrIndex=strfind(TxtFile,BriteTarget);

for i=1:length(BriteStrIndex)
    Name(i,:)=TxtFile([BriteStrIndex(i)+length(BriteTarget):BriteStrIndex(i)+length(BriteTarget)+6]); 
end
Name=table(Name);
BriteStrIndex(length(BriteStrIndex)+1)=length(TxtFile);
SampleAll=0*ones(size(Name,1),size(GeneAbundance,2));
 
 for i=1:size(GeneName,1)
     i
    GeneStrIndex=strfind(TxtFile,[char(GeneName{i,:}) '</a></dt>']);
    for j=1:length(GeneStrIndex)
        for k=1:length(BriteStrIndex)-1
            if GeneStrIndex(j)>BriteStrIndex(k)&&GeneStrIndex(j)<BriteStrIndex(k+1)
                Fitness(k,i)=log(Abun(i))-log(Phi(i))-mm;
                break;
            end
        end
    end
 end
 
 save('Fitness.mat','Fitness','Name');
 for i=1:size(Fitness,1)
    for j=1:size(Fitness,2)
        if Fitness(i,j)==-Inf||Fitness(i,j)==0
            Fitness(i,j)=NaN;    
        end
    end
 end
 tt=Fitness(2:end,:);
 nn=Name{2:end,1};
 bbmm=nanmean(tt,2);
 [B I]=sort(bbmm,'descend');
 boxplot(tt(I,:)');