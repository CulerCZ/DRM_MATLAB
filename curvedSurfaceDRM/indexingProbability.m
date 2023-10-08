%% indexing rate of a masked simulated patterns

% using indexing method on the simulated DRPs, which are basically indexing
% distance calculation
outcome_180 = cell(length(5:5:45),2);
for tiltAngle = 5:5:45
    mask = createTiltMask(tiltAngle,exp_para);
    mask_double = repmat(mask,[1,2]);
    mask_180 = mask_double(:,37:108);
    % create masked DRP dictionary
    drpDic_masked = maskDRPLib(drpLib.drpDic,mask_180);
    fprintf("Dictionary of %02d tilting masked DRPs is established!\n",tiltAngle);

    % perform indexing on masked simulated DRPs
    nn = length(drpDic_masked);
    euler = zeros(nn,3);
    K = 3;
    idx = zeros(nn,K);
    for ii = 1:nn
        [euler(ii,:),idx(ii,:)] = DirectDIEngine(drpDic_masked(ii), drpLib, K=K);
        workbar(ii/nn,sprintf("processing DRP %d / %d",[ii nn]));
    end
    outcome_180{tiltAngle/5,1} = euler;
    outcome_180{tiltAngle/5,2} = idx;
    fprintf("Indexing of %02d masking is finished!\n",tiltAngle);
end
%% test the indexing speed with respect to the different numbers of K
batchNum = nn/2;
totalIter = ceil(nn / batchNum);
tic
for ii = 1:totalIter
    idx_temp = ((ii-1)*batchNum + 1 : min(ii*batchNum,nn))';
    [euler(idx_temp,:),idx(idx_temp,:)] = DirectDIEngine(drpDic_masked(idx_temp), drpLib, K=K);
    workbar(idx_temp(end)/nn,sprintf("processing DRP %d / %d",[idx_temp(end) nn]));
end
toc
%% Choose batch number equals nn/2 as final selection
thetaRange = 0:5:60;
outcome = cell(length(thetaRange),2);
phi = repmat(linspace(exp_para.ph_min, exp_para.ph_max, exp_para.ph_num),exp_para.th_num,1);
theta = repmat(linspace(exp_para.th_min, exp_para.th_max, exp_para.th_num)',1,exp_para.ph_num);
ii = 0;
for tiltAngle = thetaRange
    % mask = createTiltMask(tiltAngle,exp_para);
    ii = ii+1;
    mask = checkVecAvailable(phi, theta, 180, tiltAngle);
    % mask_180 = checkVecAvailable(phi, theta, 0, tiltAngle);

    % create masked DRP dictionary
    drpDic_masked = maskDRPLib(drpLib.drpDic,mask);
    fprintf("Dictionary of %02d tilting masked DRPs is established!\n",tiltAngle);

    % perform indexing on masked simulated DRPs
    nn = length(drpDic_masked);
    euler = zeros(nn,3);
    K = 3;
    idx = zeros(nn,K);
    batchNum = nn/2;
    [euler(1:batchNum,:),idx(1:batchNum,:)] = DirectDIEngine( ...
        drpDic_masked(1:batchNum), drpLib, K=K);
    [euler(batchNum+1:nn,:),idx(batchNum+1:nn,:)] = DirectDIEngine( ...
        drpDic_masked(batchNum+1:nn), drpLib, K=K);
    outcome{tiltAngle/5+1,1} = euler;
    outcome{tiltAngle/5+1,2} = idx;
    fprintf("Indexing of %02d masking is finished!\n",tiltAngle);
end

%% postprocessing data

ori_ref = orientation.byEuler(drpLib.eulerDic.*degree,cs);
numPoints = size(outcome,1);
angleDiffTilt = zeros(nn,numPoints);
for ii = 1:numPoints
    ori_temp = orientation.byEuler(outcome{ii,1}.*degree,cs);
    angleDiffTilt(:,ii) = angle(ori_ref,ori_temp,cs)./degree;
end

figure("Units","centimeters","Position",[2 2 20 8])
plot(thetaRange,sum(angleDiffTilt > 5)/nn,':','Marker','o','LineWidth',2,'Color','black')
hold on
% plot(sum(angleDiffTilt > 10)/nn,'LineWidth',2)
% plot(sum(angleDiffTilt > 15)/nn,'LineWidth',2)
% hold off
set(gca,'LineWidth',2,'FontSize',14)
xticks(0:9)
xticklabels(0:5:45)
ylim([0 0.6])
yticks(0:0.1:0.6)
yticklabels(0:10:60)

plot(0:9,1-sum(misOriPixel < 15)/size(misOriPixel,1)/0.9,':','Marker','diamond','LineWidth',2,'Color','black')

peakEff = [0 2.23 4.45 6.4 9.02 11.76 14.37 16.72 19.73]./100;
bandEff = [0 25.43 36.45 42.56 46.27 49.90 51.87 54.49 57.09]./100;

plot(peakEff,'--','Marker',"v",'LineWidth',2,'Color','black')
plot(bandEff,'-.','Marker','square','LineWidth',2,'Color','black')

legend("Indexing Loss (SIM)","Indexing Loss (EXP)","Peak Loss","Band Loss",'Location','northeastoutside')
xlabel("Tilt angle (deg)",'FontSize',14)
ylabel("Information Loss (%)",'FontSize',14)
% exportgraphics(gcf,"/Users/chenyangzhu/Desktop/datasetTemp/Tilt_InfoLoss.tif",'Resolution',600)

%% plot the orientation failed in IPF figure
figure("Units","centimeters","Position",[2 2 24 24])
tiledlayout(3,3)

for ii = 1:9
    failedCases = angleDiffTilt(:,ii) > 5;
    orientationFailed = ori_ref(failedCases);
    n_failed = length(orientationFailed);
    colorData = ones(n_failed,1).*angleDiffTilt(failedCases,ii)/15;
    nexttile(ii)
    plotDataIPF(colorData,orientationFailed,vector3d.Z,"markersize",4);
    % title(sprintf("Tilt angle %d degree", ii*5),'FontSize',16)
end
% exportgraphics(gcf,"/Users/chenyangzhu/Desktop/datasetTemp/Tilt_InfoLoss_IPF.tif",'Resolution',600)

figure("Units","centimeters","Position",[2 2 24 24])
tiledlayout(3,3)
for ii = 1:9
    failedCases = angleDiffTilt(:,ii) > 5;
    orientationFailed = ori_ref(failedCases);
    n_failed = length(orientationFailed);
    % colorData = ones(n_failed,1).*angleDiffTilt(failedCases,ii)/15;
    nexttile(ii)
    histogram(angleDiffTilt(failedCases,ii),'BinLimits',[1 62],'BinWidth',2)
    set(gca,"LineWidth",2,"FontSize",14)
    % title(sprintf("Tilt angle %d degree", ii*5),'FontSize',16)
end
% exportgraphics(gcf,"/Users/chenyangzhu/Desktop/datasetTemp/Tilt_InfoLoss_histo.tif",'Resolution',300)


ori_ref = orientation.byEuler(drpLib.eulerDic.*degree,cs);
angleDiffTilt = zeros(nn,9);
setMTEXpref('xAxisDirection','east'); 
setMTEXpref('zAxisDirection','outOfPlane'); 
for ii = 3:9
    ori_temp = orientation.byEuler(outcome{ii,1}.*degree,cs);
    angleDiffTilt(:,ii) = angle(ori_ref,ori_temp,cs)./degree;
    selectCase = angleDiffTilt(:,ii) > 5;
    selectNum = sum(selectCase);
    selectOri = ori_temp(selectCase);
    selectangleDiff = angleDiffTilt(selectCase,ii);
    % odf = selectangleDiff*unimodalODF(selectOri,'halfwidth',1*degree);
    % figure, plotPDF(odf,Miller({0 0 1},cs),'noTitle','resolution',1*degree)
    % figure, plotPDF(selectOri,selectangleDiff,Miller({0 0 1},cs),'noTitle','all','MarkerSize',2)
    figure, plotPDF(ori_ref,angleDiffTilt(:,ii),Miller({0 0 1},cs),'noTitle','all','MarkerSize',max(1,10*angleDiffTilt(:,ii)/15))
    clim([0 15])
    colormap(flipud(gray))
    % colorbar
    % exportgraphics(gcf,sprintf("/Users/chenyangzhu/Desktop/datasetTemp/tilt_%02d_ErrorFunc.tif",ii*5),'Resolution',600)
    close gcf
end
%% ------
ori_ref = orientation.byEuler(drpLib.eulerDic.*degree,cs);
angleDiffTilt = zeros(nn,9);
setMTEXpref('xAxisDirection','east'); 
setMTEXpref('zAxisDirection','outOfPlane'); 
for ii = 3:9
    ori_temp = orientation.byEuler(outcome{ii,1}.*degree,cs);
    angleDiffTilt(:,ii) = angle(ori_ref,ori_temp,cs)./degree;
    selectCase = angleDiffTilt(:,ii) > 5;
    selectNum = sum(selectCase);
    selectOri = ori_temp(selectCase);
    selectangleDiff = angleDiffTilt(selectCase,ii);
    % odf = selectangleDiff*unimodalODF(selectOri,'halfwidth',1*degree);
    % figure, plotPDF(odf,Miller({0 0 1},cs),'noTitle','resolution',1*degree)
    % figure, plotPDF(selectOri,selectangleDiff,Miller({0 0 1},cs),'noTitle','all','MarkerSize',2)
    figure, plotPDF(ori_ref,angleDiffTilt(:,ii),Miller({0 0 1},cs),'noTitle','all','MarkerSize',max(1,10*angleDiffTilt(:,ii)/15))
    clim([0 15])
    colormap(flipud(gray))
    % colorbar
    exportgraphics(gcf,sprintf("/Users/chenyangzhu/Desktop/datasetTemp/tilt_%02d_ErrorFunc.tif",ii*5),'Resolution',600)
    close gcf
end



%% the fraction of data available for different reflectance features
% it will be a probability function as a function of crystallographic
% information as well as material

% the probability function should be plotted in IPF possibly

% y = linspace(exp_para.th_min,exp_para.th_max,exp_para.th_num);
% x = linspace(exp_para.ph_min,exp_para.ph_max,exp_para.ph_num);
% [phi,theta] = meshgrid(x,y);
% phi0 = 180;
% theta0 = 20;
% dirAvailable = checkVecAvailable(phi, theta, phi0, theta0);
% figure, imshow(dirAvailable)

% iterating through the DRP dictionary and calculate the fraction of
% reflectance features loss

outcome_10 = indexMaskedDRPs(drpLib,10,exp_para);
%%
performance = procOutcome(outcome_10,drpLib);

%% processing outcome and acquire probability mapping on the inverse pole figure
ori_refPoint = orientation.byEuler(euDic_R.*degree,cs);
r = inv(ori_refPoint) * vector3d.Z;
xyz_Point = normr(sort(abs(r.xyz),2,'ascend'));
xy_Point = sort(xyz_Point(:,1:2),2,'descend');
resolution = 1000;
x = linspace(0,sqrt(2)/2,resolution);
y = linspace(0,0.6,resolution);
[xx,yy] = meshgrid(x,y);
plotFlag = true;
% plotIdx = [1 2 8 20 38 61 89 120 154 190];
plotIdx = 38:2:60;
for idx = 1:length(plotIdx)
    % outcome_10{ii,3} = angle(ori_refPoint,orientation.byEuler(outcome_10{ii,1}.*degree,cs),cs)./degree;
    ii = plotIdx(idx);
    z = outcome_10{ii,3};
    if plotFlag
        zz = griddata(xy_Point(:,1),xy_Point(:,2),z,xx,yy,"linear");
        figure, surf(xx,yy,zz,'FaceColor','interp','EdgeColor','none')
        view(0,90)
        xlim([0 1])
        ylim([0 0.6])
        set(gca, 'Visible','off')
        axis equal
        clim([0 15])
        hold on
        plot(boundaryX,boundaryY,'k',LineWidth=6)
        exportgraphics(gcf,sprintf("/Users/chenyangzhu/Desktop/datasetTemp/TiltIndexErrorFig/%03d_%03d.tif",[floor(rhoGrid(ii)), floor(thetaGrid(ii))]),"Resolution",300)
        close gcf
    end
end


%%
poleValue = zeros(length(rhoGrid),3);
poleValue(:,1) = rhoGrid';
poleValue(:,2) = thetaGrid';
poleValue(:,3) = mean(performance.angDiff);
poleX = cosd(poleValue(:,2)).*cosd(poleValue(:,1));
poleY = cosd(poleValue(:,2)).*sind(poleValue(:,1));
poleZ = sind(poleValue(:,2));
poleXX = poleX ./(1+poleZ);
poleYY = poleY ./(1+poleZ);
resolution = 100;
x = linspace(-1,1,resolution);
y = linspace(-1,1,resolution);
[xx,yy] = meshgrid(x,y);
poleZZ = griddata(poleXX,poleYY,poleValue(:,3),xx,yy,'linear');

% poleZZ = griddata(poleXX,poleYY,sum(errOutcome>5)/9866,xx,yy,'linear');
figure("Units","centimeters","Position",[2 2 16 16])
surf(xx,yy,poleZZ,'FaceColor','interp','EdgeColor','none')

clim([0,15])
% clim([0 0.5])
colormap("parula")
view(0,90)
xlim([-1 1])
ylim([-1 1])
axis equal
set(gca,'Visible','off')

hold on
t = linspace(0,2*pi,360);
plot(0.98*cos(t),0.98*sin(t),'color','black','LineWidth',6)
%% Function to be used in this script
function mask = createTiltMask(thetaThres, exp_para)
% create tilted self-blockage mask to the raw data
    y = linspace(exp_para.th_min,exp_para.th_max,exp_para.th_num);
    x = linspace(exp_para.ph_min,exp_para.ph_max,exp_para.ph_num);
    [phph,thth] = meshgrid(x,y);
    theta_critical = asind(abs(sind(thetaThres)*sind(atand(-1./(tand(phph)*cosd(thetaThres))))));
    theta_critical(:,x<90|x>270) = 0;
    mask = thth > theta_critical;
end

function drpDic_masked = maskDRPLib(drpDic,mask)
    nn = length(drpDic);
    drpDic_masked = cell(nn,1);
    for ii = 1:nn
        drp_temp = drpDic{ii};
        drp_temp(~mask) = 0;
        drpDic_masked{ii} = drp_temp;
    end
end

function [Euler, idx] = DirectDIEngine(drpM, drpLib, options)
    arguments
        drpM
        drpLib
        options.K (1,1) double = 1
    end
    drplist_s = zeros(length(drpLib.drpDic), numel(drpM{1,1}));
    for ii = 1:length(drpLib.drpDic)
        drplist_s(ii,:) = double(reshape(drpLib.drpDic{ii},1,[]))/256;
    end

    n2 = length(drpM);
    drplist_m = zeros(n2, numel(drpM{1,1}));
    for jj = 1:n2
        drplist_m(jj,:) = double(reshape(drpM{jj},1,[]))/256;
    end
    % drplist_m = double(reshape(drpM{1},1,[]))/256;

    [idx, ~] = knnsearch(drplist_s, drplist_m, K=options.K);
    Euler = drpLib.eulerDic(idx(:,1),:);
end


function dirAvailable = checkVecAvailable(phi, theta, phi0, theta0, options)
    % function to check the available of faceting direction available or
    % not with knowing polar position of peak direction (phi, theta)
    arguments
        phi double
        theta double
        phi0 (1,1) double = 0
        theta0 (1,1) double = 0
        options.forbiddenRegion
    end
    nPhi = length(phi);
    if length(theta) ~= nPhi
        error("the numbers of phi and theta are not consistent.")
    end
    % dirAvailable = false(nPhi,1);
    cond1 = tand(theta)./sqrt(cosd(phi).^2*cosd(phi0)^2+sind(phi).^2*sind(phi0)^2) > tand(theta0);
    cond2 = cosd(phi)*cosd(phi0) + sind(phi)*sind(phi0) < 0;
    dirAvailable = cond1 | cond2;
end


function outcome = indexMaskedDRPs(drpLib,resolution,exp_para,options)
    arguments
        drpLib
        resolution
        exp_para
        options.K (1,1) double = 3
    end
    drpDic = drpLib.drpDic;
    S2G = equispacedS2Grid('resolution',resolution*degree,'maxTheta',90*degree);
    rhoGrid = S2G.rho ./ degree;
    thetaGrid = 90 - S2G.theta ./ degree;
    
    % create DRP angular positions
    y = linspace(exp_para.th_min,exp_para.th_max,exp_para.th_num);
    x = linspace(exp_para.ph_min,exp_para.ph_max,exp_para.ph_num);
    [phi,theta] = meshgrid(x,y);
    
    outcome = cell(length(rhoGrid),2);
    timeIter = zeros(length(rhoGrid),1);
    
    for ii = 1:length(rhoGrid)
        tic
        mask = checkVecAvailable(phi, theta, 0, thetaGrid(ii));
        shift = floor(rhoGrid(ii)/5);
        mask = circshift(mask,shift,2);
        drpDic_masked = maskDRPLib(drpDic,mask);
        % fprintf("Dictionary of (%03d / %03d) tilting masked DRPs is established!\n",[floor(rhoGrid(ii)), floor(thetaGrid(ii))]);
    
        nn = length(drpDic_masked);
        euler = zeros(nn,3);
        K = options.K;
        idx = zeros(nn,K);
        batchNum = nn/2;
        [euler(1:batchNum,:),idx(1:batchNum,:)] = DirectDIEngine( ...
            drpDic_masked(1:batchNum), drpLib, K=K);
        [euler(batchNum+1:nn,:),idx(batchNum+1:nn,:)] = DirectDIEngine( ...
            drpDic_masked(batchNum+1:nn), drpLib, K=K);
        outcome{ii,1} = euler;
        outcome{ii,2} = idx;
        fprintf("Finish %04d / %04d tilt case.\n",[ii length(rhoGrid)]);
        timeIter(ii) = toc;
        totalSecLeft = (length(rhoGrid)-ii)*mean(timeIter(1:ii));
        totalMinLeft = totalSecLeft/60;
        fprintf("Estimate time remaining: %d mins %02d seconds.\n", ...
            [floor(totalMinLeft), floor((totalMinLeft-floor(totalMinLeft))*60)])
    end

end

function performance = procOutcome(outcome,drpLib,options)
    arguments
        outcome
        drpLib
        options.crystalSymmetry = crystalSymmetry('cubic')
    end
    cs = options.crystalSymmetry;
    ori_ref = orientation.byEuler(drpLib.eulerDic.*degree,cs);
    nn = length(outcome);
    performance.angDiff = zeros(length(drpLib.drpDic),nn);
    performance.oriRef = ori_ref;
    for ii = 1:length(drpLib.drpDic)
        performance.odf(ii) = unimodalODF(ori_ref(ii),'halfwidth',3*degree);
        workbar(ii/length(drpLib.drpDic))
    end
    for ii = 1:nn
        ori_temp = orientation.byEuler(outcome{ii,1}.*degree,cs);
        performance.angDiff(:,ii) = angle(ori_ref, ori_temp, cs)./degree;
    end

end