degrees = 0:5:45;
num_deg = length(degrees);
misOriGrain_data = cell(1,num_deg);
for ii = 1:num_deg
    data_temp = load(fullfile("/Users/chenyangzhu/Desktop/datasetTemp", ...
        [num2str(degrees(ii),'%02d'),'deg'],"/indexed.mat"),"misOriAngle");
    data00 = data_temp.misOriAngle;
    misOriGrain_data{ii} = data00(~isnan(data00));
    workbar(ii/num_deg,sprintf("load %02d / %02d folders",[ii num_deg]));
end
clear data00 data_temp
misOriPixel = cell2mat(misOriGrain_data);

%%
imgFolder = "/Users/chenyangzhu/Desktop/datasetTemp";
misOri_03 = misOriPixel(misOriPixel(:,1)>30,:);
% figure, bar3(misOri_03)

misOriGrainAve = mean(misOriPixel,2);
misOriGrainStd = std(misOriPixel,1,2);
grainFeature = [misOriGrainAve, misOriGrainStd, max(misOriPixel,[],2),min(misOriPixel,[],2)];

grainTotOff = grainFeature(:,4) > 20;
grainWell = grainFeature(:,3) < 20;
grainJump = grainFeature(:,2) > 16;
grainStable = grainFeature(:,2) < 3;


figure("Units","centimeters","Position",[2 2 4*sum(grainTotOff)/10 4])
pcolor(misOriPixel(grainTotOff,:)')
clim([0 62])
axis equal
ylim([1 10])
xlim([1 sum(grainTotOff)])
xticks([])
yticks([])
colormap(turbo)
exportgraphics(gca,fullfile(imgFolder,"grainTotOff.tiff"),'Resolution',300)
%%
figure("Units","centimeters","Position",[2 2 4*sum(grainWell)/10 4])
pcolor(misOriPixel(grainWell,:)')
clim([0 20])
axis equal
ylim([1 10])
xticks([])
yticks([])
colormap(turbo)
exportgraphics(gca,fullfile(imgFolder,"grainWell.tiff"),'Resolution',300)
%%
figure("Units","centimeters","Position",[2 2 4*sum(grainJump)/10 4])
pcolor(misOriPixel(grainJump,:)')
clim([0 62])
axis equal
ylim([1 10])
xticks(1:sum(grainJump))
xticklabels(find(grainJump))
% xtickangle(60)
yticks([])
colormap(jet)
exportgraphics(gca,fullfile(imgFolder,"grainJump.tiff"),'Resolution',300)

%%
figure("Units","centimeters","Position",[2 2 4*sum(grainStable)/10 4])
pcolor(misOriPixel(grainStable,:)')
clim([0 62])
axis equal
ylim([1 10])
xticks([])
yticks([])
colormap(turbo)
exportgraphics(gca,fullfile(imgFolder,"grainStable.tiff"),'Resolution',300)

%% ranking by standard deviation
% [~,Idx] = sort(grainFeature(:,2));
Idx = 1:size(grainFeature,1);
for ii = 1:5
    observation = misOriPixel(Idx((ii-1)*90+1:ii*90),:);
    figure("Units","centimeters","Position",[2 2 20 4])
    pcolor(observation')
    clim([0 62])
    axis equal
    ylim([1 10])
    yticks([])
    % xticks(linspace(1,90,10))
    % xticklabels((0:10:90) + (ii-1)*90)
    % xtickangle(45)
    xticks([])
    colormap(jet)
    exportgraphics(gca,fullfile(imgFolder,sprintf("grainCat_%01d.tiff",ii)),'Resolution',300)
end
%% grain category mapping
IdxTotOff = Idx_whole(grainTotOff);
IdxWell = Idx_whole(grainWell);
IdxgrainJump = Idx_whole(grainJump);
grainCat = zeros(size(grainIdMap_ebsd));
grainCat(grainIdMap_ebsd>0) = 1;
for ii = 1:length(IdxTotOff)
    grainCat(grainIdMap_ebsd == IdxTotOff(ii)) = 2;
end
for ii = 1:length(IdxWell)
    grainCat(grainIdMap_ebsd == IdxWell(ii)) = 3;
end
for ii = 1:length(IdxgrainJump)
    grainCat(grainIdMap_ebsd == IdxgrainJump(ii)) = 4;
end
figure, imshow(grainCat,[0 4.2],"Border","tight")
cm = [228,26,28;55,126,184;77,175,74;152,78,163;255,127,0]./255;
colormap(cm)
% colorbar("FontSize",14)
exportgraphics(gcf,fullfile(imgFolder,"grainCategory.tiff"),'Resolution',300)
%%
figure, plotIPDF(ori_mean(IdxTotOff),repmat(cm(2,:),length(IdxTotOff),1),vector3d.Z,"MarkerSize",80*grains.grainSize(IdxTotOff))
hold on
plotIPDF(ori_mean(IdxWell),repmat(cm(3,:),length(IdxWell),1),vector3d.Z,"MarkerSize",80*grains.grainSize(IdxWell))
plotIPDF(ori_mean(IdxgrainJump),repmat(cm(4,:),length(IdxgrainJump),1),vector3d.Z,"MarkerSize",80*grains.grainSize(IdxgrainJump))
exportgraphics(gcf,fullfile(imgFolder,"grainCatIPFZ.tiff"),'Resolution',300)

%%
figure, hold on
for ii = 1:size(misOriGrainMat,1)
    if misOriGrainMat(ii,1) > 15 && misOriGrainMat(ii,10) > 15
        plot(degrees, misOriGrainMat(ii,:))
    else
        continue
    end
end

%%
a = prctile(misOriPixel,25);
b = prctile(misOriPixel,50);
c = prctile(misOriPixel,75);
figure('unit','centimeters','position',[1 1 20 15])
yyaxis left
% plot(x,a,'LineWidth',2,'Color','#fdae61')
hold on
plot(degrees,b,'LineWidth',2,'Color','#3182bd')
% plot(x,c,'LineWidth',2,'Color','#fdae61')
fill([degrees, fliplr(degrees)], [a, fliplr(c)], [158,202,225]/255, 'FaceAlpha', 0.3, 'EdgeColor', [158,202,225]/255);
xlabel("tilting angle (deg)")
ylabel("indexing error (deg)")
box on
set(gca,'LineWidth',2,'FontSize',14)

idxRate = [87.34 84.28 82.57 80.98 77.21 68.50 62.00 53.75 44.94 38.49];
yyaxis right
plot(degrees, idxRate, 'LineWidth', 2, 'Marker', '+')
ylabel("indexing rate (%)", "FontSize", 14)
