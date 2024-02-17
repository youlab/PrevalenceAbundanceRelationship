clear;
clc;
Nt=100;
m=1;
d=10^[-5];
CC=linspecer(2);
p=0.08;

xs=0:0.01:0.4;
for hjk=1:length(xs)
    x=xs(hjk);
    p1(hjk)=x^(Nt*m*p-1)*(1-x)^(Nt*m*(1-p)-1);
    c=gamma(Nt)/gamma(Nt*(1-p))/gamma(Nt*p);
    p1(hjk)=c*p1(hjk)*0.01;
    p2(hjk)=nchoosek(Nt,fix(x*Nt))*p^(x*Nt)*(1-p)^(Nt-Nt*x);
    
end
subplot(1,2,1);
plot(xs,p1,'-','markersize',10,'color',CC(1,:),'linewidth',2);hold on;
plot(xs,p2,'-','markersize',10,'color',CC(2,:),'linewidth',2);hold on;
set(gca,'fontsize',12);
xlabel('abundance','fontsize',16);
ylabel('probability density','fontsize',16);

subplot(1,2,2);
plot(p1,p2,'ko','markersize',10,'linewidth',1);hold on;
set(gca,'fontsize',12);
xlabel('prob beta','fontsize',16);
ylabel('prob binomial','fontsize',16);

set(gcf,'position',[100 100 600 270]);

saveas(gcf,'ProbCompare.fig');
saveas(gcf,'ProbCompare.pdf');