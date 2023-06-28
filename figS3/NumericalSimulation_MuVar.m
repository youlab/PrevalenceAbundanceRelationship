clear;
clc;
global NumSpecies D sigma gamma mu Nm;
NumSpecies=100;
Nm=ones(NumSpecies,1);
gamma=-ones(NumSpecies,NumSpecies);    
sigma=0*ones(NumSpecies,1);

D=0.02;

NumWell=100;
NumCell=100;
NumRound=10;
GrowthTime=0:24;

temp=rand(NumSpecies,1);
InitialProb=temp/sum(temp);

varrs=[0 0.05 0.1 0.15];
C=linspecer(length(varrs));
Slope=0;

for jfk=1:length(varrs)
    Prob=InitialProb;
    varr=varrs(jfk);
    mu=0.1+varr*rand(NumSpecies,1);

    PreAbund=0*ones(NumSpecies,NumWell,NumRound);
    PostAbund=0*ones(NumSpecies,NumWell,NumRound);

    for round=1:NumRound  
        round+(jfk-1)*NumRound
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
        ff=0.8;
        Prob=ff*Prob+(1-ff)*temp(1:NumSpecies)/sum(temp(1:NumSpecies));
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
        Slope(jfk,i,j)=p(1);
        end
    end
    
end

for i=1:length(varrs)
    plot([1:NumRound],Slope(i,:,2),'o-','color',C(i,:),'markersize',8,'linewidth',1.5);hold on;
end
axis([1 NumRound 0.9 max(Slope(:))+0.05]);
set(gca,'fontsize',16);
xlabel('cycles','fontsize',20);
ylabel('scaling coefficient \theta','fontsize',20);
legend({'','','','',''});
legend boxoff;
set(gcf,'position',[100 100 250 250]);
saveas(gcf,'NumericalSimulation_MuVar.fig');
save('NumericalSimulation_MuVar.mat');

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
