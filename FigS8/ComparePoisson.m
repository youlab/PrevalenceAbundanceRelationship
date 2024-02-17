clear;
clc;
close all;
ms=[5,10,20,50];
ra=0:0.01:1;
CC=linspecer(length(ms));

gamma=0.1;

for j=1:length(ms)
    m=ms(j);
    for i=1:length(ra)
        a=fix(m*ra(i));
        p1(i)=nchoosek(m,a)*gamma^a*(1-gamma)^(m-a);
        p2(i)=exp(-m*gamma)*(m*gamma)^a/factorial(a);
    end
    plot(p1,p2,'-','linewidth',3,'color',CC(j,:));hold on;
end
plot(p1,p1,'k--','linewidth',0.5);
axis([10^(-5),1,10^(-5),1]);
set(gca,'XScale','log');
set(gca,'YScale','log');
xticks(10.^[-5:0]);
yticks(10.^[-5:0]);
set(gca,'fontsize',8);
xlabel('probability in random sampling');
ylabel('probability predicted by Poisson distribution');

set(gcf,'position',[100,100,250,250]);

saveas(gcf,'ComparePoisson.fig');
saveas(gcf,'ComparePoisson.pdf');