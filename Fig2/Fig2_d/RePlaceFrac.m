clear;
clc;
close all;
C=linspecer(2);
NumGene=100;
NumWell=100;
NumCell=200;
NumRound=20;
GrowthTime=24;
Mu=0+0.1*rand(NumGene,1);
temp=rand(NumGene,1);
InitialProb=temp/sum(temp);

ffs=[0.05 0.1 0.2 0.5];
C=linspecer(length(ffs));
Slope=0;
for jfk=1:length(ffs)
    Prob=InitialProb;
    ff=ffs(jfk);
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
        Slope(jfk,i,j)=p(1);
        end
    end
end


subplot(1,2,1);
for i=1:length(ffs)
    plot([1:2:NumRound],Slope(i,1:2:end,2),'o-','color',C(i,:),'markersize',5,'linewidth',1.5);hold on;
end
axis([1 NumRound 0.9 max(Slope(:))+0.05]);
set(gca,'fontsize',16);
xlabel('cycles','fontsize',20);
ylabel('scaling coefficient \theta','fontsize',20);
legend({'','','','',''});
legend boxoff;
ylim([0.9 1.4]);
set(gcf,'position',[100 100 500 250]);
saveas(gcf,'ReplaceFrac.fig');
saveas(gcf,'ReplaceFrac.pdf');

