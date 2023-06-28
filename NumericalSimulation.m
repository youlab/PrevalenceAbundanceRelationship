clear;
clc;
global NumSpecies D sigma gamma mu Nm;
NumSpecies=100;
Nm=ones(NumSpecies,1);
gamma=-ones(NumSpecies,NumSpecies);    
sigma=0*ones(NumSpecies,1);
mu=0.1+0.1*rand(NumSpecies,1);

D=0.02;

NumWell=100;
NumCell=100;
NumRound=4;
GrowthTime=0:24;

temp=rand(NumSpecies,1);
InitialProb=temp/sum(temp);
Prob=InitialProb;

PreAbund(NumSpecies,NumWell,NumRound)=0;
PostAbund(NumSpecies,NumWell,NumRound)=0;

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
save('NumericalSimulation.mat');
NumericalSimulation_Visual
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
