clear;
clc;
close all;
C=linspecer(2);
NumGene=100;
NumWells=[20 50 100 200];
NumCell=200;
NumRound=3;
GrowthTime=24;
Mu=0+0.1*rand(NumGene,1);


RR=20;
Slope=0;
ff=1;
for iop=1:length(NumWells)
    NumWell=NumWells(iop);
    for fgh=1:RR
        fgh
        temp=rand(NumGene,1);
        InitialProb=temp/sum(temp);
        Prob=InitialProb;
        PreAbund=0*ones(NumGene,NumWell,NumRound);
        PostAbund=0*ones(NumGene,NumWell,NumRound);
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
            Prob=ff*temp/sum(temp)+(1-ff)*Prob;
        end
    
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
                abund(k)=mean((Data(k,:)));
            end
            phi=-log(1-prev);
    
            A=log(phi((prev<1)&(prev>0)));
            B=log(abund((prev<1)&(prev>0)));
            p=polyfit(A,B,1);
            Slope(i,j,fgh,iop)=p(1);
            end
        end
    end
end


CC=linspecer(4);
for i=1:length(NumWells)
plot(NumWells(i)-3+6*rand(RR,1),reshape(Slope(1,1,:,i),[RR 1]),'o','markersize',5,'color',0.7*[1 1 1]);hold on;
end
h1=errorbar(NumWells,reshape(mean(Slope(1,1,:,:),3),[length(NumWells),1]),reshape(std(Slope(1,1,:,:),0,3),[length(NumWells),1]),'color',CC(1,:),'linewidth',3);
hold on;

for i=1:length(NumWells)
plot(NumWells(i)-3+6*rand(RR,1),reshape(Slope(3,2,:,i),[RR 1]),'o','markersize',5,'color',0.7*[1 1 1]);hold on;
end
h2=errorbar(NumWells,reshape(mean(Slope(3,2,:,:),3),[length(NumWells),1]),reshape(std(Slope(3,2,:,:),0,3),[length(NumWells),1]),'color',CC(2,:),'linewidth',3);
hold on;
xlim([0 max(NumWells)+10]);

% legend([h1,h2],{'',''});
% legend boxoff;
% [h,p]=ttest(reshape(Slope(1,1,:,1),[RR 1]),reshape(Slope(3,2,:,1),[RR 1]))
set(gca,'fontsize',16);
xlabel('Number of communities','fontsize',20);
ylabel('\theta value','fontsize',20);
set(gcf,'position',[100 100 270 270]);
saveas(gcf,'ThetaVariations.fig');
saveas(gcf,'ThetaVariations.pdf');

