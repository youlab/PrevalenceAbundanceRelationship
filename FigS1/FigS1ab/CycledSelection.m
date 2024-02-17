clear;
clc;
close all;
C=linspecer(2);
NumGene=100;
NumWell=100;
NumCell=100;
NumRound=5;
GrowthTime=24;
Mu=0.1*rand(NumGene,1);
temp=rand(NumGene,1);
InitialProb=temp/sum(temp);
Prob=InitialProb;

PreAbund(NumGene,NumWell,NumRound)=0;
PostAbund(NumGene,NumWell,NumRound)=0;
for round=1:NumRound
    round
    tape=0;
    tape(length(Prob)+1)=0;
    for i=2:length(tape)
        tape(i)=tape(i-1)+Prob(i-1);
    end
    for i=1:NumWell
        
        tt=NumCell;%poissrnd(NumCell);
        for j=1:tt
            dd=rand;
            for k=1:length(tape)-1
                if dd>tape(k)&&dd<=tape(k+1)
                    PreAbund(k,i,round)=PreAbund(k,i,round)+1;
                    break;
                end
            end
        end  
    end
    
    for i=1:NumGene
        for j=1:NumWell
            PostAbund(i,j,round)=PreAbund(i,j,round)*exp(Mu(i)*GrowthTime);
        end
    end
    
    temp=sum(PostAbund(:,:,round),2);
    temp=temp/sum(temp);
    temp=temp*1+Prob*0;
    Prob=temp/sum(temp);
end

CCR=zeros(NumRound,1);

for i=1:NumRound
    for j=1
        if j==1
            Data=PreAbund(:,:,i);
        end
        if j==2
            Data=PostAbund(:,:,i);
        end
        Data=Data./sum(Data,1);
    
        prev=0;
        abund=0;
        phi=0;
        for k=1:size(Data,1)
            prev(k)=nnz(Data(k,:))/size(Data,2);
        end
        
        CCR(i)=corr(Mu,prev','type','Spearman');

    end

    figure(1);
    subplot(2,3,i);
    plot(Mu,prev,'o','Color',C(1,:),'Markersize',4);hold on;
    set(gca,'fontsize',8);

end

figure(1);
subplot(2,3,4);
xlabel('Replication rate','fontsize',12);
ylabel('Occupancy','FontSize',12);
set(gcf,'position',[100,100,600,400]);
saveas(gcf,'CycledSelection_1.fig');
saveas(gcf,'CycledSelection_1.pdf');


figure(2);
plot(1:NumRound,CCR,'o-','markersize',10);
set(gcf,'position',[100,100,220,220]);
set(gca,'fontsize',8);
xlabel('cycles','fontsize',20);
ylabel('Spearman correlation','FontSize',20);
axis([0.8 5 -0.05 1]);
saveas(gcf,'CycledSelection_2.fig');
saveas(gcf,'CycledSelection_2.pdf');


