% exp_para_111 = exp_para;
% exp_para_111.faceting = [1 1 1];
% exp_para_111.fitting_para = [1 0.6 15 4 0.8 8];
% drpLib_111 = DRPLibGenerator(3*degree,exp_para_111);
thetaRange = 0:1:50;
outcome_111_masked = cell(length(thetaRange),2);
phi = repmat(linspace(exp_para.ph_min, exp_para.ph_max, exp_para.ph_num),exp_para.th_num,1);
theta = repmat(linspace(exp_para.th_min, exp_para.th_max, exp_para.th_num)',1,exp_para.ph_num);
ii = 0;
for tiltAngle = thetaRange
    % mask = createTiltMask(tiltAngle,exp_para);
    ii = ii+1;
    mask = createTiltMask(tiltAngle,exp_para);
    % create masked DRP dictionary
    drpDic_masked = maskDRPLib(drpLib_111.drpDic,mask);
    drpLib_111_masked = drpLib_111;
    drpLib_111_masked.drpDic = maskDRPLib(drpLib_111.drpDic,mask);
    fprintf("Dictionary of %02d tilting masked DRPs is established!\n",tiltAngle);

    % perform indexing on masked simulated DRPs
    nn = length(drpDic_masked);
    euler = zeros(nn,3);
    K = 3;
    idx = zeros(nn,K);
    batchNum = nn/2;
    % drpLib_111_masked = maskDRPLib(drpLib_111)
    [euler(1:batchNum,:),idx(1:batchNum,:)] = DirectDIEngine( ...
        drpDic_masked(1:batchNum), drpLib_111_masked, K=K);
    [euler(batchNum+1:nn,:),idx(batchNum+1:nn,:)] = DirectDIEngine( ...
        drpDic_masked(batchNum+1:nn), drpLib_111_masked, K=K);
    outcome_111_masked{tiltAngle+1,1} = euler;
    outcome_111_masked{tiltAngle+1,2} = idx;
    fprintf("Indexing of %02d masking is finished!\n",tiltAngle);
end

%% plot the curve
nn = 22892;
% ori_ref = orientation.byEuler(drpLib.eulerDic.*degree,cs);
numPoints = size(outcome,1);
angleDiffTilt = zeros(nn,numPoints);
angleDiffTilt_111 = angleDiffTilt;
for ii = 1:numPoints
    ori_temp = orientation.byEuler(outcome_111{ii,1}.*degree,cs);
    angleDiffTilt(:,ii) = angle(ori_ref,ori_temp,cs)./degree;
    ori_temp = orientation.byEuler(outcome_111_masked{ii,1}.*degree,cs);
    angleDiffTilt_111(:,ii) = angle(ori_ref,ori_temp,cs)./degree;
end
clear ori_temp

figure("Units","centimeters","Position",[2 2 18 9])
plot(1:51,sum(angleDiffTilt > 5)/nn,'-','LineWidth',2,'Color','red')
hold on
plot(1:51,sum(angleDiffTilt_111>5)/nn,'-','LineWidth',2,'Color','blue')
% plot(0:5:45,1-sum(misOriPixel < 20)/size(misOriPixel,1)/0.9, ...
%     ':','Marker','o','LineWidth',2,'Color','black')

% legend("Sim_{(100)}","Sim_{(111)}","Exp",'Location','northwest')
legend("Sim_{(111)-shadowDI}","Sim_{(111)}","Exp",'Location','northwest')

xlabel("Tilt angle (deg)",'FontSize',14)
ylabel("Indexing Loss (%)",'FontSize',14)

set(gca,'LineWidth',2,'FontSize',14)
xticks(0:5:45)
xticklabels(0:5:45)
xlim([0 45])
ylim([0 0.6])
yticks(0:0.1:0.6)
yticklabels(0:10:60)
exportgraphics(gcf,fullfile(saveFolder,"indexingSimPerformance_shadowDI.tif"),'Resolution',600)


%% need an function to identify the shadowed region
drp_temp = drpDic_masked{1000};
figure, imshow(drp_temp < 10)
shadowRegion = findShadowRegion(drp_temp, 10, boundary=1);
figure, imshow(shadowRegion)

%% algorithm to accelerate the indexing speed




%% functions used in this script
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


function shadowRegion = findShadowRegion(drp, thres, options)
    arguments
        drp
        thres (1,1) uint8 = 10
        options.boundary (1,1) double = 1
    end
    % [n1,n2] = size(drp);
    shadowRegion = drp < thres;
    if options.boundary > 0
        perimeter = bwperim(shadowRegion);
        se = strel('square',3);
        dilatedPerimeter = imdilate(perimeter,se);
        shadowRegion(dilatedPerimeter & shadowRegion) = 0;
    end
    
end


