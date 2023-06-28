clear;
clc;
close all;
load('StrainAbundances.mat');
load('GrowthRate.mat');
C=linspecer(2);
    
for i=1:4
    if i==1
        Data=C0;
    end
    if i==2
        Data=C1;
    end
    if i==3
        Data=C2;
    end
    if i==4
        Data=C3;
    end
    
    Data=Data./sum(Data,1);
    for k=1:length(Data(:))
        if isnan(Data(k))==1
            Data(k)=0;
        end
    end
    Prev=0;
    Phi=0;
    Abund=0;
    for k=1:size(Data,1)
        Prev(k)=nnz(Data(k,:))/size(Data,2);
        Phi(k)=-log(1-Prev(k));
        Abund(k)=mean((Data(k,:)));            
    end
    figure(i);
    if i==1
        plot(log(Phi),log(Abund),'.','color',[132,212,212]/256,'markersize',30);hold on;
    else
        plot(log(Phi),log(Abund),'.','color',C(2,:),'markersize',30);hold on;
    end
    x=log(Phi);
    y=log(Abund);
    
    A=x(Prev>0&Prev<1);
    B=y(Prev>0&Prev<1);
    p=polyfit(A,B,1);
    p
    xxx=1;
    if i==1
        plot([min(A)-xxx max(A)+xxx],p(1)*[min(A)-xxx max(A)+xxx]+p(2),'-','color',[132,212,212]/256,'linewidth',3);hold on;
    else
        plot([min(A)-xxx max(A)+xxx],p(1)*[min(A)-xxx max(A)+xxx]+p(2),'-','color',C(2,:),'linewidth',3);hold on;
        plot([min(A)-xxx max(A)+xxx],[min(A)-xxx max(A)+xxx]+mean(B-A),'--','color',[132,212,212]/256,'linewidth',3);hold on;
    end
    set(gca,'fontsize',12);
    xlabel('log\phi','fontsize',16);
    ylabel('log\alpha','fontsize',16); 
    if i==1
        axis([-6 3 -12 -3]);
    end
    if i==2
        axis([-6 3 -12 -3]);
    end
    if i==3
        axis([-6 3 -15 -4]);
    end
    if i==4
        axis([-6 3 -15 -3]);
    end
    set(gcf,'position',[100 100 260 260]);
end


