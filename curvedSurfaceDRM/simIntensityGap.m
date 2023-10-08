%% the gap between simulation data and experimental data
figure,

scatter(0:5:45,resultGap,72,'black','filled')
%%
figure,
yyaxis left
scatter(0:5:45,resultGap,72,"blue",'filled')
yyaxis right
scatter(0:5:45,sind(45+(0:5:45))/sind(45),72,"red",'filled')

%%
resultSim = [0,sum(angleDiffTilt>5)/nn];
resultExp = 1-sum(misOriPixel<15)/size(misOriPixel,1)/0.9;
resultGap = (resultExp-resultSim)*100;
% P = polyfit(0:5:45, resultGap, 4);
% yfit = polyval(P,0:45);
kfit = (sind(45:5:90)-sind(45))'\resultGap';
figure("Units","centimeters","Position",[2 2 20 8])
scatter(sind(45:5:90),resultGap,72,'black','filled')
hold on
box on
xfit = linspace(sind(45),sind(90),50);
plot(xfit, kfit*(xfit-sind(45)), 'k-.','LineWidth',2)
set(gca,'LineWidth',2,'FontSize',14)

xticks(sind(45:5:90))
ylim([0 30])
yticks(0:5:30)
legend("Discrepency","linear fitting",'Location','northeastoutside')
xlabel("Tilt angle \alpha (deg)",'FontSize',14)
ylabel("Indexing Loss (%)",'FontSize',14)
% exportgraphics(gcf,"/Users/chenyangzhu/Desktop/datasetTemp/ExpSimGap.tif",'Resolution',600)
