clear;
clc;
close all;
load('NumericalSimulation.mat');
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

Data=PreAbund(:,:,1);
Data=Data./sum(Data,1);
prev=0;
abund=0;
phi=0;
for k=1:size(Data,1)
    prev(k)=nnz(Data(k,:))/size(Data,2);
    abund(k)=mean((Data(k,:)));
end
phi=-log(1-prev);
figure(1);
plot(log(phi),log(abund),'.','markersize',20,'color',[132,212,212]/256);hold on;

A=log(phi((prev<1)&(prev>0)));
B=log(abund((prev<1)&(prev>0)));
p=polyfit(A,B,1);
Slope=p(1);
plot([min(A)-0.2 max(A)+0.2],p(1)*[min(A)-0.2 max(A)+0.2]+p(2),'--','color',[132,212,212]/256,'linewidth',3);hold on;

set(gca,'fontsize',12);
xlabel('log\phi','fontsize',16);
ylabel('log\alpha','fontsize',16);
axis([-4.5 2 -10 -4]);

set(gcf,'position',[100 100 260 260]);
saveas(gcf,sprintf('NumericalSimulationVisual_0.fig'));
saveas(gcf,sprintf('NumericalSimulationVisual_0.eps'));
 
Slope=0;
for i=1:NumRound
    for j=[1 2]
        if j==1
            Data=PreAbundGene(:,:,i);
        end
        if j==2
            Data=PostAbundGene(:,:,i);
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
    figure(i+1);
    if j==2
            plot(log(phi),log(abund),'.','markersize',20,'color',[236,84,84]/256);hold on;
    end
    A=log(phi((prev<1)&(prev>0)));
    B=log(abund((prev<1)&(prev>0)));
    p=polyfit(A,B,1);
    Slope(i,j)=p(1);
    if j==2
    plot([min(A)-0.2 max(A)+0.2],p(1)*[min(A)-0.2 max(A)+0.2]+p(2),'--','color',C(2,:),'linewidth',3);hold on;
    plot([min(A)-0.5 max(A)+0.5],[min(A)-0.5 max(A)+0.5]+mean(B-A),'--','color',[132,212,212]/256,'linewidth',3);hold on;
    end
%     text(0.1,0.9, sprintf('slope = %f', p(1)),'fontsize',16,'Units','normalized','Color','b');
    set(gca,'fontsize',12);
        if j==1
            xlabel('log\phi','fontsize',16);
            ylabel('log\alpha','fontsize',16);
        end
    end
    if i==1
        axis([-4.5 2 -10.5 -1.5]);
    end
    if i==2
        axis([-3 2 -11 -3]);
    end
    if i==3
        axis([-3 2 -12 -4]);
    end
    set(gcf,'position',[100 100 260 260]);
    saveas(gcf,sprintf('NumericalSimulationVisual_%.0f.fig',i));
    saveas(gcf,sprintf('NumericalSimulationVisual_%.0f.eps',i));
end