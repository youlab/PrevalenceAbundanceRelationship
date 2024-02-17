clear;
clc;
close all;
C=linspecer(2);
NumGene=100;
NumWell=100;
NumCell=100;
NumRound=100;
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
    temp=temp*0.03+Prob*0.97;
    Prob=temp/sum(temp);
end

Slope=zeros(NumRound,2);
Chi=zeros(NumRound,NumGene);
GenePhi=zeros(NumRound,NumGene);
GeneAbun=zeros(NumRound,NumGene);
for i=1:NumRound
    for j=[1 2]
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
            abund(k)=mean(Data(k,:));
        end
        phi=-log(1-prev);
        A=log(phi((prev<1)&(prev>0)));
        B=log(abund((prev<1)&(prev>0)));
        p=polyfit(A,B,1);
        Slope(i,j)=p(1);
        if j==2
            Chi(i,:)=log(abund)-log(phi)-mean(B-A);
            GenePhi(i,:)=phi;
            GeneAbun(i,:)=abund;
        end
        % subplot(10,NumRound/10,i);
        % plot(A,B,'k.');hold on;
    end
end

figure(1);
plot(1:NumRound,Slope(:,2),'o','markersize',5);
set(gcf,'position',[100,100,220,220]);
set(gca,'fontsize',8);
xlabel('cycles','fontsize',20);
ylabel('\theta','FontSize',20);
xlim([0 100]);
saveas(gcf,'CycledSelection_1.fig');
saveas(gcf,'CycledSelection_1.pdf');


figure(2);
boxplot(Chi','BoxStyle','filled','Colors',C(1,:));
set(gca,'fontsize',8);
set(gcf,'position',[100 100 1200 220]);
xlabel('cycles','fontsize',20);
ylabel('\chi','fontsize',20);
saveas(gcf,'CycledSelection_2.fig');
saveas(gcf,'CycledSelection_2.pdf');

figure(3);
for i=1:NumRound
    temp=Chi(i,:);
    MM=nanmean(temp(temp~=-Inf));
    SS=nanstd(temp(temp~=-Inf));
    errorbar(Slope(i,2),MM,SS,'o','markersize',5,'capsize',0,'color',C(1,:));hold on;
end
set(gca,'fontsize',8);
set(gcf,'position',[100 100 620 220]);
xlabel('\theta','fontsize',20);
ylabel('\chi','fontsize',20);
saveas(gcf,'CycledSelection_3.fig');
saveas(gcf,'CycledSelection_3.pdf');