clear;
clc;
global NumSpecies D sigma gamma mu Nm;
NumSpecies=100;
Nm=ones(NumSpecies,1);
gamma=-ones(NumSpecies,NumSpecies);    
sigma=0*ones(NumSpecies,1);

mu=0.1+0.1*rand(NumSpecies,1);

D=0.02;

NumWells=[20 50 100 200];
RR=20;

NumCell=100;
NumRound=3;
GrowthTime=0:24;

Slope=0*ones(NumRound,2,RR,length(NumWells));

for iop=1:length(NumWells)
    NumWell=NumWells(iop);
    for fgh=1:RR
        fgh
        temp=rand(NumSpecies,1);
        InitialProb=temp/sum(temp);
        Prob=InitialProb;

        PreAbund=0*ones(NumSpecies,NumWell,NumRound);
        PostAbund=0*ones(NumSpecies,NumWell,NumRound);

        for round=1:NumRound  
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
                initial=PreAbund(:,i,round)/NumCell/10^6;
                [t,y]=ode45(@multi_species,GrowthTime,initial);

                for k=1:NumSpecies
                    PostAbund(k,i,round)=y(end,k);
                end
            end

            temp=sum(PostAbund(:,:,round),2);
            
            ff=0.1/NumSpecies;
            Prob=ff*ones(NumSpecies,1)+temp(1:NumSpecies)/sum(temp(1:NumSpecies));
            Prob=Prob/sum(Prob);

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
save('NumericalSimulation.mat');

CC=linspecer(4);
for i=1:length(NumWells)
plot(NumWells(i)-3+6*rand(RR,1),reshape(Slope(1,1,:,i),[RR 1]),'o','markersize',5,'color',0.7*[1 1 1]);hold on;
end
h1=errorbar(NumWells,reshape(mean(Slope(1,1,:,:),3),[length(NumWells),1]),reshape(std(Slope(1,1,:,:),0,3),[length(NumWells),1]),'color',CC(1,:),'linewidth',3);
hold on;
% errorbar(NumWells,reshape(mean(Slope(1,2,:,:),3),[length(NumWells),1]),reshape(std(Slope(1,2,:,:),0,3),[length(NumWells),1]),'color',CC(2,:),'linewidth',2);
% hold on;
% errorbar(NumWells,reshape(mean(Slope(2,2,:,:),3),[length(NumWells),1]),reshape(std(Slope(2,2,:,:),0,3),[length(NumWells),1]),'color',CC(3,:),'linewidth',2);
% hold on;
for i=1:length(NumWells)
plot(NumWells(i)-3+6*rand(RR,1),reshape(Slope(3,2,:,i),[RR 1]),'o','markersize',5,'color',0.7*[1 1 1]);hold on;
end
h2=errorbar(NumWells,reshape(mean(Slope(3,2,:,:),3),[length(NumWells),1]),reshape(std(Slope(3,2,:,:),0,3),[length(NumWells),1]),'color',CC(2,:),'linewidth',3);
hold on;
xlim([10 210]);
% legend([h1,h2],{'Neutral','Deterministic'});
legend([h1,h2],{'',''});
legend boxoff;
% [h,p]=ttest(reshape(Slope(1,1,:,1),[RR 1]),reshape(Slope(3,2,:,1),[RR 1]))
set(gca,'fontsize',16);
xlabel('Number of communities','fontsize',20);
ylabel('\theta value','fontsize',20);
set(gcf,'position',[100 100 270 270]);
saveas(gcf,'NumericalSimulation.fig');
saveas(gcf,'NumericalSimulation.eps');


function dydt=multi_species(t,y)
        global NumSpecies D sigma gamma mu Nm;
        dydt(NumSpecies,1)=0;
        for i=1:NumSpecies
            Neg=0;
            Pos=0;
            for j=1:NumSpecies
                if gamma(j,i)<0
                Neg=Neg-gamma(j,i)*y(j);
                end
                if gamma(j,i)>0
                Pos=Pos+gamma(j,i)*y(j);
                end
            end
            
            mui=1-Neg/Nm(i)+sigma(i)*Pos/(Pos+1)/Nm(i);
     
            dydt(i,1)=mu(i)*y(i)*mui-D*y(i);
        end      
end
