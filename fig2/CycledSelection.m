clear;
clc;
close all;
C=linspecer(2);
NumGene=100;
NumWell=100;
NumCell=100;
NumRound=3;
GrowthTime=24;
Mu=0.1*rand(NumGene,1);
temp=rand(NumGene,1);
InitialProb=temp/sum(temp);
Prob=InitialProb;

PreAbund(NumGene,NumWell,NumRound)=0;
PostAbund(NumGene,NumWell,NumRound)=0;
for round=1:NumRound
    tape=0;
    tape(length(Prob)+1)=0;
    for i=2:length(tape)
        tape(i)=tape(i-1)+Prob(i-1);
    end
    for i=1:NumWell
        i
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
    Prob=temp/sum(temp);
end

Slope=0;
Data=PreAbund(:,:,1);
prev=0;
abund=0;
phi=0;
for k=1:size(Data,1)
    prev(k)=nnz(Data(k,:))/size(Data,2);
    abund(k)=mean((Data(k,:)));
end
phi=-log(1-prev);
plot(log(phi),log(abund),'.','markersize',30,'color',[132,212,212]/256);hold on;
A=log(phi((prev<1)&(prev>0)));
B=log(abund((prev<1)&(prev>0)));
p=polyfit(A,B,1);
plot([min(A)-0.2 max(A)+0.2],p(1)*[min(A)-0.2 max(A)+0.2]+p(2),'-','color',[132,212,212]/256,'linewidth',3);hold on;
set(gca,'fontsize',12);
xlabel('log\phi','fontsize',16);
ylabel('log\alpha','fontsize',16);
axis([-4.1 2 -4.1 1.5]);
set(gcf,'position',[100 100 260 260]);
saveas(gcf,'CycledSelection_1.fig');
saveas(gcf,'CycledSelection_1.pdf');

Slope=0;
Residue=0;
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
        figure(i+1);
        A=log(phi((prev<1)&(prev>0)));
        B=log(abund((prev<1)&(prev>0)));
        p=polyfit(A,B,1);
        Slope(i,j)=p(1);
        if j==2
            plot(log(phi),log(abund),'.','markersize',30,'color',C(2,:));hold on;
            plot([min(A)-2 max(A)+2],p(1)*[min(A)-2 max(A)+2]+p(2),'-','color',C(2,:),'linewidth',3);hold on;
            plot([min(A)-2 max(A)+2],[min(A)-2 max(A)+2]+mean(B-A),'--','color',[132,212,212]/256,'linewidth',3);hold on;
        end
        
        set(gca,'fontsize',12);
        xlabel('log\phi','fontsize',16);
        ylabel('log\alpha','fontsize',16);
        set(gcf,'position',[100 100 260 260]);
        axis([-5 2 -12 -2]);
        saveas(gcf,sprintf('CycledSelection_%.0f.fig',i+1));
        saveas(gcf,sprintf('CycledSelection_%.0f.pdf',i+1));
    end
end

figure(i+2);
Data=PostAbund(:,:,2);
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
plot(Mu,log(abund)-log(phi)-mean(B-A),'.','markersize',30,'color',C(2,:));hold on;
axis([-0.01 0.11 -1.3 1.3]);
set(gca,'fontsize',14);
xlabel('\mu','fontsize',16);
ylabel('log\alpha-log\phi','fontsize',16);
set(gcf,'position',[100 100 260 260]);
saveas(gcf,sprintf('CycledSelection_%.0f.fig',i+2));
saveas(gcf,sprintf('CycledSelection_%.0f.pdf',i+2));

save('CycledSelection.mat');
Slope