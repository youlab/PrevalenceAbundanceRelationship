clear;
clc;
close all;
load('Fitness.mat');
C=linspecer(3);
for i=1:size(Fitness,1)
    for j=1:size(Fitness,2)
        if Fitness(i,j)==-Inf||Fitness(i,j)==0
            Fitness(i,j)=NaN;    
        end
    end
 end
 tt=Fitness(2:end,:);
 nn=Name{2:end,1};
 bbmm=nanmean(tt,2);
 [B I]=sort(bbmm);
 boxplot(tt(I,:)','symbol','','colors',[255 87 51]/256,'orientation', 'horizontal');hold on;
 axis([-3 3 0.5 53.5]);
 plot([0 0],[0.5 53.5],'--','linewidth',1,'color',C(3,:));hold on;
 set(gcf,'position',[100 100 200 500]);
 nnl=nn(I(end:-1:1),:);
 saveas(gcf,'FitnessVisualization.fig');
 saveas(gcf,'FitnessVisualization.png');