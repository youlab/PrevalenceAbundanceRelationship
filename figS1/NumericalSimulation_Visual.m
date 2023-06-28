clear;
clc;
close all;
load('NumericalSimulation_PosFreq01.mat');
C=linspecer(2);

for i=NumRound
    for j=2
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
    figure(i+1);
    if j==2
            plot(log(phi),log(abund),'.','markersize',30,'color',C(2,:));hold on;
    end
    A=log(phi((prev<1)&(prev>0)));
    B=log(abund((prev<1)&(prev>0)));
    p=polyfit(A,B,1);
    Slope(i,j)=p(1);
    if j==2
        plot([min(A)-0.5 max(A)+0.5],p(1)*[min(A)-0.5 max(A)+0.5]+p(2),'-','color',C(2,:),'linewidth',3);hold on;
        plot([min(A)-0.5 max(A)+0.5],[min(A)-0.5 max(A)+0.5]+mean(B-A),'--','color',[132,212,212]/256,'linewidth',3);hold on;
    end
%     text(0.1,0.9, sprintf('slope = %f', p(1)),'fontsize',16,'Units','normalized','Color','b');
    set(gca,'fontsize',12);
        if j==1
            xlabel('log\phi','fontsize',16);
            ylabel('log\alpha','fontsize',16);
        end
    end
    axis([-4.5 2.5 -10.5 -1.5]);
    set(gcf,'position',[100 100 260 260]);
    saveas(gcf,'NumericalSimulationVisual_PosFreq01.fig');
    saveas(gcf,'NumericalSimulationVisual_PosFreq01.png');
end

Slope