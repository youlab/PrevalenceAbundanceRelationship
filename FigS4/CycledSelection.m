clear;
clc;
close all;
C=linspecer(2);
NumGenes=10:10:200;
NumCell=100;
NumWell=100;
NumRepeat=500;
Slope=NaN*ones(NumRepeat,length(NumGenes));
for dfg=1:length(NumGenes)
    dfg
    NumGene=NumGenes(dfg);
    PreAbund=zeros(NumGene,NumWell,NumRepeat);
    for repeat=1:NumRepeat
        Prob=rand(NumGene,1);
        Prob=Prob/sum(Prob);
        for i=1:NumWell
            for j=1:NumCell
                dd=rand;
                
                tape=zeros(NumGene+1,1);
                for ii=2:length(tape)
                    tape(ii)=tape(ii-1)+Prob(ii-1);
                end
    
                for k=1:length(tape)-1
                    if dd>tape(k)&&dd<=tape(k+1)
                        PreAbund(k,i,repeat)=PreAbund(k,i,repeat)+1;
                        break;
                    end
                end
            end  
        end
        
    end
    
    for as=1:NumRepeat
        Data=PreAbund(:,:,as);
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
        p=polyfit(A,B,1);
        if length(A)>=3 && abs(p(1))<10
            Slope(as,dfg)=p(1);
        end
    end
end
figure(1);
errorbar(NumGenes,nanmean(Slope,1),nanstd(Slope,1),'o','markersize',6,'color',C(1,:),'linewidth',2,'CapSize',0);
xlim([7,max(NumGenes)+5]);
set(gcf,'position',[100 100 250 250]);
set(gca,'fontsize',12);
xlabel('Number of genes','fontsize',14);
ylabel('\theta values','fontsize',14);
saveas(gcf,'CycledSelection_1.fig');
saveas(gcf,'CycledSelection_1.pdf');

figure(2);
plot(NumGenes,nanstd(Slope,1),'o-','markersize',6,'color',C(1,:),'linewidth',2);hold on;
plot(1:200,0.02*ones(200,1),'--','markersize',6,'color',C(1,:),'linewidth',1);
xlim([7,max(NumGenes)+5]);
set(gcf,'position',[100 100 250 250]);
set(gca,'fontsize',12);
set(gca,'YScale','log');
xlabel('Number of genes','fontsize',14);
ylabel('\theta std','fontsize',14);
saveas(gcf,'CycledSelection_2.fig');
saveas(gcf,'CycledSelection_2.pdf');