clear;
close all;
clc;
load('MyColormap.mat');
load('Ocean\RawDataGene.mat');
Data=TotalReads(2:end,:);
Data=Data./sum(Data,1);
Prev=0;
Abun=0;
for j=1:size(Data,1)
    Prev(j)=nnz(Data(j,:))./size(Data,2);
    Abun(j)=mean((Data(j,:)));
end
figure(1);
Phi=-log(1-Prev);
x=log(Phi(Prev>0&Prev<1));
y=log(Abun(Prev>0&Prev<1));
p=polyfit(x,y,1)
q=corr(x',y','type','Pearson')
XXX=[x',y'];
N=hist3(XXX,'CDataMode','auto','FaceColor','interp','Ctrs',{-5:0.05:2 -24:0.1:-8});
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(-5,2,size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(-24,-8,size(N_pcolor,1)); % Rows of N_pcolor
h = pcolor(xl,yl,N_pcolor);
set(h, 'EdgeColor', 'none');
c=mymap;
colormap(c);
caxis([0 5]);
hold on;
plot([min(x) 1.2*max(x)],p(1)*[min(x) 1.2*max(x)]+p(2),'-','color','k','linewidth',2.5);hold on;
plot([min(x) 1.2*max(x)],[min(x) 1.2*max(x)]+mean(y-x),'--','color','k','linewidth',2.5);hold on;
set(gca,'fontsize',16);
xlabel('log\phi','fontsize',20);
ylabel('log\alpha','fontsize',20);
set(gcf,'position',[100 100 350 350]);
saveas(gcf,'fig4_Genes_Ocean.fig');
saveas(gcf,'fig4_Genes_Ocean.png');


load('HumanSkin\RawDataGeneNose.mat');
Data=SampleHumanSkin;
Data=Data./sum(Data,1);
Prev=0;
Abun=0;
for j=1:size(Data,1)
    Prev(j)=nnz(Data(j,:))./size(Data,2);
    Abun(j)=mean((Data(j,:)));
end
figure(2);
Phi=-log(1-Prev);
x=log(Phi(Prev>0&Prev<1));
y=log(Abun(Prev>0&Prev<1));
p=polyfit(x,y,1)
q=corr(x',y','type','Pearson')
XXX=[x',y'];
N=hist3(XXX,'CDataMode','auto','FaceColor','interp','Ctrs',{-5:0.05:2 -24:0.1:-6});
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(-5,2,size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(-24,-6,size(N_pcolor,1)); % Rows of N_pcolor
h = pcolor(xl,yl,N_pcolor);
set(h, 'EdgeColor', 'none');
c=mymap;
colormap(c);
caxis([0 7]);
hold on;
plot([min(x) 1.2*max(x)],p(1)*[min(x) 1.2*max(x)]+p(2),'-','color','k','linewidth',2.5);hold on;
plot([min(x) 1.2*max(x)],[min(x) 1.2*max(x)]+mean(y-x),'--','color','k','linewidth',2.5);hold on;
set(gca,'fontsize',16);
xlabel('log\phi','fontsize',20);
ylabel('log\alpha','fontsize',20);
set(gcf,'position',[100 100 350 350]);
saveas(gcf,'fig4_Genes_HumanNose.fig');
saveas(gcf,'fig4_Genes_HumanNose.png');

load('HumanSkin\RawDataGeneCheek.mat');
Data=SampleHumanSkin;
Data=Data./sum(Data,1);
Prev=0;
Abun=0;
for j=1:size(Data,1)
    Prev(j)=nnz(Data(j,:))./size(Data,2);
    Abun(j)=mean((Data(j,:)));
end
figure(3);
Phi=-log(1-Prev);
x=log(Phi(Prev>0&Prev<1));
y=log(Abun(Prev>0&Prev<1));
p=polyfit(x,y,1)
q=corr(x',y','type','Pearson')
XXX=[x',y'];
N=hist3(XXX,'CDataMode','auto','FaceColor','interp','Ctrs',{-5:0.05:2 -24:0.1:-6});
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(-5,2,size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(-24,-6,size(N_pcolor,1)); % Rows of N_pcolor
h = pcolor(xl,yl,N_pcolor);
set(h, 'EdgeColor', 'none');
c=mymap;
colormap(c);
caxis([0 7]);
hold on;
plot([min(x) 1.2*max(x)],p(1)*[min(x) 1.2*max(x)]+p(2),'-','color','k','linewidth',2.5);hold on;
plot([min(x) 1.2*max(x)],[min(x) 1.2*max(x)]+mean(y-x),'--','color','k','linewidth',2.5);hold on;
set(gca,'fontsize',16);
xlabel('log\phi','fontsize',20);
ylabel('log\alpha','fontsize',20);
set(gcf,'position',[100 100 350 350]);
saveas(gcf,'fig4_Genes_HumanCheek.fig');
saveas(gcf,'fig4_Genes_HumanCheek.png');

load('HumanGut\RawDataGene.mat');
Data=SampleHuman;
Data=Data./sum(Data,1);
Prev=0;
Abun=0;
for j=1:size(Data,1)
    Prev(j)=nnz(Data(j,:))./size(Data,2);
    Abun(j)=mean((Data(j,:)));
end
figure(4);
Phi=-log(1-Prev);
x=log(Phi(Prev>0&Prev<1));
y=log(Abun(Prev>0&Prev<1));
p=polyfit(x,y,1)
q=corr(x',y','type','Pearson')
XXX=[x',y'];
N=hist3(XXX,'CDataMode','auto','FaceColor','interp','Ctrs',{-5:0.05:2 -21:0.1:-9});
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(-5,2,size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(-21,-9,size(N_pcolor,1)); % Rows of N_pcolor
h = pcolor(xl,yl,N_pcolor);
set(h, 'EdgeColor', 'none');
c=mymap;
colormap(c);
caxis([0 7]);
hold on;
plot([min(x) 1.2*max(x)],p(1)*[min(x) 1.2*max(x)]+p(2),'-','color','k','linewidth',2.5);hold on;
plot([min(x) 1.2*max(x)],[min(x) 1.2*max(x)]+mean(y-x),'--','color','k','linewidth',2.5);hold on;
set(gca,'fontsize',16);
xlabel('log\phi','fontsize',20);
ylabel('log\alpha','fontsize',20);
set(gcf,'position',[100 100 350 350]);
saveas(gcf,'fig4_Genes_HumanGut.fig');
saveas(gcf,'fig4_Genes_HumanGut.png');

load('Pig\RawDataGene.mat');
Data=SamplePig;
Data=Data./sum(Data,1);
Prev=0;
Abun=0;
for j=1:size(Data,1)
    Prev(j)=nnz(Data(j,:))./size(Data,2);
    Abun(j)=mean((Data(j,:)));
end
figure(5);
Phi=-log(1-Prev);
x=log(Phi(Prev>0&Prev<1));
y=log(Abun(Prev>0&Prev<1));
p=polyfit(x,y,1)
q=corr(x',y','type','Pearson')
XXX=[x',y'];
N=hist3(XXX,'CDataMode','auto','FaceColor','interp','Ctrs',{-5:0.05:2 -20:0.1:-8});
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(-5,2,size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(-20,-8,size(N_pcolor,1)); % Rows of N_pcolor
h = pcolor(xl,yl,N_pcolor);
set(h, 'EdgeColor', 'none');
c=mymap;
colormap(c);
caxis([0 5]);
hold on;
plot([min(x) 1.2*max(x)],p(1)*[min(x) 1.2*max(x)]+p(2),'-','color','k','linewidth',2.5);hold on;
plot([min(x) 1.2*max(x)],[min(x) 1.2*max(x)]+mean(y-x),'--','color','k','linewidth',2.5);hold on;
set(gca,'fontsize',16);
xlabel('log\phi','fontsize',20);
ylabel('log\alpha','fontsize',20);
set(gcf,'position',[100 100 350 350]);
saveas(gcf,'fig4_Genes_PigGut.fig');
saveas(gcf,'fig4_Genes_PigGut.png');

load('Mouse\RawDataGene.mat');
Data=SampleMouse;
Data=Data./sum(Data,1);
Prev=0;
Abun=0;
for j=1:size(Data,1)
    Prev(j)=nnz(Data(j,:))./size(Data,2);
    Abun(j)=mean((Data(j,:)));
end
figure(6);
Phi=-log(1-Prev);
x=log(Phi(Prev>0&Prev<1));
y=log(Abun(Prev>0&Prev<1));
p=polyfit(x,y,1)
q=corr(x',y','type','Pearson')
XXX=[x',y'];
N=hist3(XXX,'CDataMode','auto','FaceColor','interp','Ctrs',{-5:0.05:2 -24:0.1:-6});
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(-5,2,size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(-24,-6,size(N_pcolor,1)); % Rows of N_pcolor
h = pcolor(xl,yl,N_pcolor);
set(h, 'EdgeColor', 'none');
c=mymap;
colormap(c);
caxis([0 4]);
hold on;
plot([min(x) 1.2*max(x)],p(1)*[min(x) 1.2*max(x)]+p(2),'-','color','k','linewidth',2.5);hold on;
plot([min(x) 1.2*max(x)],[min(x) 1.2*max(x)]+mean(y-x),'--','color','k','linewidth',2.5);hold on;
set(gca,'fontsize',16);
xlabel('log\phi','fontsize',20);
ylabel('log\alpha','fontsize',20);
set(gcf,'position',[100 100 350 350]);
saveas(gcf,'fig4_Genes_MouseGut.fig');
saveas(gcf,'fig4_Genes_MouseGut.png');
