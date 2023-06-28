clear;
clc;
close all;
load('StrainAbundances.mat');
load('GrowthRate.mat');
C=linspecer(2);

   
for i=1:3
    if i==1
        Data=C1;
    end
    if i==2
        Data=C2;
    end
    if i==3
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
    temp=Prev;
    index=temp>0&temp<1;
    MM=mean(log(Abund(index))-log(Phi(index)));
    GR=mean(GrowthRate,1);
    AA=GR(index);
    BB=log(Abund(index))-log(Phi(index))-MM;
    fg=1:72;
    fitresult = fit(AA',BB','poly1');
    xx=0.5:0.0005:0.8;
    p22 = predint(fitresult,xx,0.95,'functional','on');
    curve1=transpose(p22(:,1));
    curve2=transpose(p22(:,2));
    x2 = [xx, fliplr(xx)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween,[255 128 0]/256,'FaceAlpha',0.5,'EdgeColor','none'); hold on;
        
    plot(AA,BB,'.','markersize',15,'color',C(2,:));hold on;
    set(gca,'fontsize',16);
    
    [RRR1 PPP1]=corr(AA',BB','type','Pearson');set(gca,'fontsize',16);
    Rho(i)=RRR1;
    PP(i)=PPP1;
    Pears(i)=RRR1;
    P=polyfit(AA,BB,1);
    xx=0.5:0.01:0.8;
    plot(xx,polyval(P,xx),'k--');hold on;
    if i==1
        axis([0.5 0.75 -1.2 0.8]);
    end
    if i==2
        axis([0.5 0.75 -1.2 1.1]);
    end
    if i==3
        axis([0.5 0.8 -2.5 2]);
    end
    xlabel('growth rate','fontsize',16);
    ylabel('\chi_i','fontsize',16);
    set(gcf,'position',[100 100 220 120]);
end


