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
NumRound=10;
GrowthTime=0:24;

temp=rand(NumSpecies,1);
InitialProb=temp/sum(temp);


ffs=[0];
C=linspecer(length(ffs));
Slope=0;

for jfk=1:length(ffs)
    ff=ffs(jfk);
    Prob=InitialProb;
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
        
        Prob=0.1/NumSpecies*ones(NumSpecies,1)+temp(1:NumSpecies)/sum(temp(1:NumSpecies));
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
            
            Data1=Data./sum(Data,1);

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
        
        prev1=0;
        abund1=0;
        phi1=0;
        for k=1:size(Data1,1)
            prev1(k)=nnz(Data1(k,:))/size(Data1,2);
            abund1(k)=mean((Data1(k,:)));
        end
        phi1=-log(1-prev1);

        A1=log(phi1((prev1<1)&(prev1>0)));
        B1=log(abund1((prev1<1)&(prev1>0)));
        p=polyfit(A1,B1,1);
        Slope1(jfk,i,j)=p(1);
        
        end
    end
    
end

plot(Slope(:),Slope1(:),'ko','markersize',10);hold on;
set(gca,'fontsize',16);
xlabel('without compositional effects','fontsize',20);
ylabel('with compositional effects','fontsize',20);
set(gcf,'position',[100 100 200 200]);
set(gcf,'position',[100 100 250 250]);
axis([0.9 1.7 0.9 1.7]);
plot([0.9 1.7],[0.9 1.7],'r--');hold on;
saveas(gcf,'CompositionalEffects.fig');
save('CompositionalEffects.mat');

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
