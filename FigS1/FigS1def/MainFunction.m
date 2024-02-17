clear;
clc;

CC=linspecer(4);

Nt=20;
ds=10.^[-4 -5 -6];
m=0.5;
ps=10^(-2)*[0.1:0.1:2]';
for j=1:length(ds)
    occu=0;
    d=ds(j);
    for i=1:length(ps)
        p=ps(i);
    fun = @(x) x.^(Nt*m*p-1).*(1-x).^(Nt*m*(1-p)-1);
    q = integral(fun,0,1);
    c=1/q;
    occu(i)=c*integral(fun,d,1);
    end
    plot(log(-log(1-occu(occu>0&occu<1)))',log(ps(occu>0&occu<1)),'.','markersize',20,'color',CC(j,:));hold on;
    ss=polyfit(log(-log(1-occu(occu>0&occu<1)))',log(ps(occu>0&occu<1)),1);
    slope(j)=ss(1);
end
set(gca,'fontsize',16);
% legend({'','',''});
% legend boxoff;
xlabel('log\phi','fontsize',24);
ylabel('log\alpha','fontsize',24);
% axis([0 4 -4.1 -2.2]);
set(gcf,'position',[100 100 270 270]);

saveas(gcf,'MainFunction.fig');
saveas(gcf,'MainFunction.eps');