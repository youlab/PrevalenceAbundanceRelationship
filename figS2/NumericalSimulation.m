clear;
clc;
global NumSpecies NumPlasmid Nm eta kappa D lambda sigma gamma mu;
NumSpecies=100;
NumPlasmid=NumSpecies;
NonMob=1:NumSpecies-20;

Nm=ones(NumSpecies,1);
gamma=(0.05-0.1*rand(NumSpecies,NumSpecies)).*(rand(NumSpecies,NumSpecies)>0.8);    
sigma=0.05*rand(NumSpecies,1);
mu=0.1+0.1*rand(NumSpecies,1);
for k=1:NumSpecies
       gamma(k,k)=-1;
end
        
kappamax=0.05;
kappa=0*ones(NumSpecies,NumPlasmid);
for i=1:NumPlasmid
kappa(:,i)=kappamax*rand*ones(NumSpecies,1);
end
kappa(:,NonMob)=0;
lambda=0.05-0.1*rand(NumSpecies,NumPlasmid);
lambda(:,NonMob)=0;

D=0.02;
eta=0*ones(NumPlasmid,NumSpecies,NumSpecies);
for i=1:NumPlasmid
eta(i,:,:)=0.2*rand*ones(1,NumSpecies,NumSpecies);
end
eta(NonMob,:,:)=0;

NumWell=100;
NumCell=100;
NumRound=3;
GrowthTime=0:24;

temp=rand(NumSpecies,1);
InitialProb=temp/sum(temp);
Prob=InitialProb;

PlasDistri=0*ones(NumSpecies,NumPlasmid);
for i=1:NumPlasmid
    PlasDistri(i,i)=1;
end

PreAbund(NumSpecies*(NumPlasmid+1),NumWell,NumRound)=0;
PostAbund(NumSpecies*(NumPlasmid+1),NumWell,NumRound)=0;

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
                    for kk=1:NumPlasmid
                        if rand<PlasDistri(k,kk)
                        PreAbund(NumSpecies+(k-1)*NumPlasmid+kk,i,round)=PreAbund(NumSpecies+(k-1)*NumPlasmid+kk,i,round)+1;
                        end
                    end
                    break;
                end
            end
        end    
        initial=PreAbund(:,i,round)/NumCell/10;
        [t,y]=ode45(@multi_plasmid,GrowthTime,initial);

        for k=1:NumSpecies*(NumPlasmid+1)
            PostAbund(k,i,round)=y(end,k);
        end
    end
    
    temp=sum(PostAbund(:,:,round),2);
    Prob=0.1/NumSpecies*ones(NumSpecies,1)+temp(1:NumSpecies)/sum(temp(1:NumSpecies));
    Prob=Prob/sum(Prob);
    
    for k=1:NumSpecies
        for kk=1:NumPlasmid
               PlasDistri(k,kk)=temp(NumSpecies+(k-1)*NumPlasmid+kk)/temp(k);
        end
    end
    
end

C=linspecer(2);
DetectionLimit=0;
PreAbundGene0=0*ones(NumPlasmid,NumWell,NumRound);
PostAbundGene0=0*ones(NumPlasmid,NumWell,NumRound);
PreAbundGene=0*ones(NumPlasmid,NumWell,NumRound);
PostAbundGene=0*ones(NumPlasmid,NumWell,NumRound);
for i=1:NumPlasmid
    PreAbundGene0(i,:,:)=sum(PreAbund(NumSpecies+i:NumPlasmid:end,:,:),1)./sum(PreAbund(1:NumSpecies,:,:),1);
    PreAbundGene(i,:,:)=PreAbundGene0(i,:,:).*(PreAbundGene0(i,:,:)>DetectionLimit);
    PostAbundGene0(i,:,:)=sum(PostAbund(NumSpecies+i:NumPlasmid:end,:,:),1)./sum(PostAbund(1:NumSpecies,:,:),1);
    PostAbundGene(i,:,:)=PostAbundGene0(i,:,:).*(PostAbundGene0(i,:,:)>DetectionLimit);
end
Slope=0;

for i=1:NumRound
    for j=[2 1]
        if j==1
            Data=PreAbundGene(:,:,i);
        end
        if j==2
            Data=PostAbundGene(:,:,i);
        end
%         Data=Data./sum(Data,1);
    
    prev=0;
    abund=0;
    phi=0;
    for k=1:size(Data,1)
        prev(k)=nnz(Data(k,:))/size(Data,2);
        abund(k)=mean((Data(k,:)));
    end
    phi=-log(1-prev);
    figure(i);
    if j==1
            plot(log(phi),log(abund),'.','markersize',20,'color',[132,212,212]/256);hold on;
    else
            plot(log(phi),log(abund),'.','markersize',20,'color',[236,84,84]/256);hold on;
    end
    A=log(phi((prev<1)&(prev>0)));
    B=log(abund((prev<1)&(prev>0)));
    p=polyfit(A,B,1);
    Slope(i,j)=p(1);
    if j==1
            plot([-0.1 1.6],p(1)*[-0.1 1.6]+p(2),'-','color',[132,212,212]/256,'linewidth',2);hold on;
    else
            plot([-0.1 1.6],p(1)*[-0.1 1.6]+p(2),'-','color',[236,84,84]/256,'linewidth',2);hold on;
    end
%     text(0.1,0.9, sprintf('slope = %f', p(1)),'fontsize',16,'Units','normalized','Color','b');
    set(gca,'fontsize',12);
        if j==1
            xlabel('log\phi','fontsize',16);
            ylabel('log\alpha','fontsize',16);
        end
    end
    set(gcf,'position',[100 100 260 260]);
    saveas(gcf,sprintf('NumericalSimulation_%.0f.fig',i));
    saveas(gcf,sprintf('NumericalSimulation_%.0f.eps',i));
end
save('NumericalSimulation.mat');

function dydt=multi_plasmid(t,y)
        global NumSpecies NumPlasmid eta kappa D lambda sigma gamma mu Nm;
        dydt(NumSpecies*(1+NumPlasmid),1)=0;
        for i=1:NumSpecies
            sum=0;
            for j=1:NumPlasmid
            sum=sum+lambda(i,j)*y(NumSpecies+NumPlasmid*(i-1)+j);
            end
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
            
            mui=mu(i)*(1-Neg/Nm(i)+sigma(i)*Pos/(Pos+1)/Nm(i));
            if y(i)==0
                dydt(i,1)=0;
            else
                dydt(i,1)=y(i)/(sum+y(i))*mui*y(i)-D*y(i);
            end
            
            for j=1:NumPlasmid
                if y(i)==0
                    betaij=0;
                elseif sum==0
                    betaij=1;
                else
                    betaij=(1+lambda(i,j))/(1+lambda(i,j)+(sum-lambda(i,j)*y(NumSpecies+NumPlasmid*(i-1)+j))/y(i));
                end
                muij=mui/(1+lambda(i,j));
                summ=0;
                for k=1:NumSpecies
                    summ=summ+eta(j,k,i)*y(NumSpecies+(k-1)*NumPlasmid+j);
                end
                dydt(NumSpecies+(i-1)*NumPlasmid+j,1)=betaij*muij*y(NumSpecies+(i-1)*NumPlasmid+j)+(y(i)-y(NumSpecies+(i-1)*NumPlasmid+j))*summ-(kappa(i,j)+D)*y(NumSpecies+(i-1)*NumPlasmid+j);
            end
        end      
end
