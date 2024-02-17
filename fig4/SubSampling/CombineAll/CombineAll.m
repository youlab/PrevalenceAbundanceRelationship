Ocean;
HumanGut;
HumanNose;
HumanCheek;
PigGut;
MouseGut;
MonkeyGut;

set(gca,'fontsize',12);
set(gcf,'position',[100 100 420 420]);
xlim([5,290]);
xlabel('sample number','fontsize',20);
ylabel('\theta value','fontsize',20);
%legend({'Ocean','Human gut','Human Nose','Human cheek','Pig gut','Mouse gut','Monkey gut'});
legend({' ',' ',' ',' ',' ',' ',' '});
legend box off;
ylim([0.9,2.4]);
saveas(gcf,'CombineAll.fig');
saveas(gcf,'CombineAll.pdf');