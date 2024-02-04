%% setting saving directory
saveFolder = "/Users/chenyangzhu/Desktop/curveDRMfigures";

%% turbine blade IPFZ
xRange = [35 605];
EUmap_select = index_result.EUmap(:,xRange(1):xRange(2),:);
figure,
imshow(plot_ipf_map(EUmap_select,plotDir="z"),'Border','tight')
exportgraphics(gcf,fullfile(saveFolder,"turbineblade_ipfz.tif"),"Resolution",300)

%% turbine blade IPFX
figure,
imshow(plot_ipf_map(EUmap_select,plotDir="x"),'Border','tight')
exportgraphics(gcf,fullfile(saveFolder,"turbineblade_ipfx.tif"),"Resolution",300)

%% turbine blade indexing error
close all
figure,
imshow(index_err_Map(:,xRange(1):xRange(2)),[0 15],'border','tight')
colormap("sky")
exportgraphics(gcf,fullfile(saveFolder,"turbineblade_idxErr.tif"),'Resolution',300)

%% plot SKY colormap colorbar
close all
plot_colorbar(cmap="sky",size=[750 50])
% exportgraphics(gcf,fullfile(saveFolder,"sky_colorbar.tif"),'Resolution',300)

%% plot horizontal line analysis of turbine blade indexing results
nn1 = size(EUmap_select,2);
y_idx = 100;
orientationList = orientation.rand(1,nn1,cs);
for ii = 1:nn1
    eaTemp = squeeze(EUmap_select(y_idx,ii,:));
    orientationList(ii) = orientation.byEuler([eaTemp(1),eaTemp(2),eaTemp(3)].*degree,cs);
end
c = jet(length(orientationList));

oriDiff = zeros(length(orientationList)-1,1);
oriCumDiff = zeros(length(orientationList)-1,1);
for ii = 1:length(oriDiff)
    oriDiff(ii) = angle(orientationList(ii+1),orientationList(ii),cs)./degree;
    oriCumDiff(ii) = angle(orientationList(ii+1),orientationList(1),cs)/degree;
end
scalingCoeff = 9.8214 / 305;  % real length / num of pixels
blockPixels = floor(1 / scalingCoeff);
numBlocks = floor(length(orientationList)/blockPixels);
numRes = length(orientationList) - blockPixels * numBlocks;

ori_pre = mean(orientationList(1:numRes));
angle_diff_list = zeros(numBlocks,4);

for ii = 1:numBlocks
    idx_temp = ((ii-1)*numBlocks+1 : ii*numBlocks)+numRes;
    ori_temp = orientationList(idx_temp);
    ori_mean = mean(ori_temp);
    ori_diff_temp = angle(ori_temp, ori_pre, cs)./degree;
    ori_diff_ori = angle(ori_temp,mean(orientationList),cs)./degree;
    angle_diff_list(ii,1) = mean(ori_diff_temp);
    angle_diff_list(ii,2) = std(ori_diff_temp);
    angle_diff_list(ii,3) = mean(ori_diff_ori);
    angle_diff_list(ii,4) = std(ori_diff_ori);
    ori_pre = ori_mean;
end
errLim = [0 10];
figure("Units",'centimeters','Position',[2 2 10 5])
errorbar(0:numBlocks-1, angle_diff_list(:,3), angle_diff_list(:,4),'Color','k','LineWidth',2)
set(gca,'LineWidth',2,'FontSize',14)
set(gca,'YLim',errLim,'XLim',[0 17])
ylabel("$\omega_{cum} (^{\circ})$",'Interpreter','latex','FontName','Monospaced','FontSize',14,'FontWeight','bold')
xlabel("x position (mm)",'FontSize',14,'FontWeight','bold')
print(gcf,fullfile(saveFolder,"topGrain_lineAnalysis"),'-dtiffn')
%% turbine blade IPDF
glabel = grainLabel(:,xRange(1):xRange(2));

oriEulerList = squeeze(reshape(EUmap_select,[],1,3));
BW_top_list = reshape(glabel,[],1)==1;
BW_bot_list = reshape(glabel,[],1)==2;
index_err_list = reshape(index_err_Map(:,xRange(1):xRange(2)),[],1);

pixel_top = orientation.byEuler(oriEulerList(BW_top_list,:).*degree,cs);
pixel_bot = orientation.byEuler(oriEulerList(BW_bot_list,:).*degree,cs);
val_top = index_err_list(BW_top_list);
val_bot = index_err_list(BW_bot_list);

odf_top = exp(-val_top) .* unimodalODF(pixel_top,'halfwidth',2*degree);
odf_bot = exp(-val_bot) .* unimodalODF(pixel_bot,'halfwidth',2*degree);
%
figure, 
plotIPDF(odf_top,zvector,'noTitle')
% hold on
% plotIPDF(orientation.byEuler(eulerangle_bot.*degree,cs),vector3d.Z,'MarkerSize',75,'MarkerFaceColor','red')
set(gcf,'Units','centimeters','Position',[2 2 8 8])
colormap("sky")
exportgraphics(gcf,fullfile(saveFolder,"topGrain_IPDF.tif"),'Resolution',300)


%% --------------------------------------------------------------------------
%% schematics of flat and tilted surface DRPs with sample
exp_para_plot.th_max = 75;
exp_para_plot.th_min = 0;
exp_para_plot.th_num = 76;
exp_para_plot.ph_min = 0;
exp_para_plot.ph_max = 359;
exp_para_plot.ph_num = 360;
exp_para_plot.faceting = [1 0 0];
exp_para_plot.fitting_para = [1, 0.6, 24, 6, 0.8, 8];
zsurf = -0.3;
V1=[-1;  1; 1; -1; -1;  1; 1; -1;];
V2=[-1; -1; 1;  1; -1; -1; 1;  1;];
V3=[-1; -1;-1; -1;  1;  1; 1;  1;]*0.1;
F= [1 2 3 4; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 5 6 7 8;];
r = 0.3;
[THETA, PHI, ~]=cart2sph(V1,V2,V3);
R=r.*ones(size(V1(:,1)));
[V1,V2,V3]=sph2cart(THETA,PHI,R);
V=[V1 V2 V3];
% apply the rotation
tilt_rot = [0 0 0];
roll = tilt_rot(1)*degree;
pitch = tilt_rot(2)*degree;
yaw = tilt_rot(3)*degree;
dcm = angle2dcm(roll, pitch, yaw, 'ZXZ');
Vertices = V*dcm;
figure("Units","centimeters","Position",[2 2 18 10])
tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
nexttile(1)
patch('Faces',F,'Vertices',Vertices,'FaceColor',[77 77 77]/255,'FaceAlpha',0.6,'EdgeColor','k',...
    'LineWidth',1); 
axis equal; hold on; view(3); % axis off;
set(gca,'visible','off')

% DRP and masked DRP
drpsim_tmp = DRPsim(0,30,45,exp_para_plot); % the matrix in form of th_num * ph_num
bgDRP_00 = createBG(0,exp_para_plot);
C = [drpsim_tmp,drpsim_tmp(:,1)];
theta = repmat(linspace(exp_para_plot.th_min, exp_para_plot.th_max, exp_para_plot.th_num)',1,exp_para_plot.ph_num+1);
phi = repmat(0:360/exp_para_plot.ph_num:360,exp_para_plot.th_num,1);
x = cosd(theta).*cosd(phi);
y = cosd(theta).*sind(phi);
z = sind(theta);

surf(x,y,zsurf*ones(size(x)),C,'EdgeColor','none','FaceAlpha',0.8)
surf(x,y,z,uint8([bgDRP_00, bgDRP_00(:,1)]/10*256),'EdgeColor','none','FaceAlpha',0.3);

axis equal
set(gca,'visible','off')
colormap('jet')
quiver3([0 0 0],[0 0 0],[0 0 0],[1.3 0 0],[0 1.3 0],[0 0 1.3],'k-','LineWidth',2)
zlim([zsurf 1.3])
xlim([-1 1.3])
ylim([-1 1.3])
view(10,10)
% tilted surfaces
nexttile(2)
bgDRP_40 = createBG(40,exp_para_plot);
bgDRP_00 = createBG(0,exp_para_plot);
bgDRP = bgDRP_00 ./ bgDRP_40;
orientation_tilt = [0,40,0];
roll = orientation_tilt(1)*degree;
pitch = orientation_tilt(2)*degree;
yaw = orientation_tilt(3)*degree;
dcm = angle2dcm(roll, pitch, yaw, 'ZYZ');
V=[V1 V2 V3];
V = V*dcm;
patch('Faces',F,'Vertices',V,'FaceColor',[77 77 77]/255,'FaceAlpha',0.6,'EdgeColor','k',...
    'LineWidth',1); axis equal; hold on; view(3);
set(gca,'visible','off')
drpsim_tmp = DRPsim(0,30,45,exp_para_plot); % the matrix in form of th_num * ph_num
drpsim_tmp(~bgDRP_40) = 0;
C = [drpsim_tmp,drpsim_tmp(:,1)];
theta = repmat(linspace(exp_para_plot.th_min, exp_para_plot.th_max, exp_para_plot.th_num)',1,exp_para_plot.ph_num+1);
phi = repmat(0:360/exp_para_plot.ph_num:360,exp_para_plot.th_num,1);
x = cosd(theta).*cosd(phi);
y = cosd(theta).*sind(phi);
z = sind(theta);

bgDRP_R = uint8([bgDRP,bgDRP(:,1)]*128);
drp_2 = uint8(double(max(C,bgDRP_R)).*double(bgDRP_R>0));
surf(x,y,zsurf*ones(size(x)),drp_2,'EdgeColor','none','FaceAlpha',0.8)
% surf(x,y,z,uint8([bgDRP_40, bgDRP_40(:,1)]/10*256),'EdgeColor','none','FaceAlpha',0.3);
surf(x,y,z,drp_2,'EdgeColor','none','FaceAlpha',0.3);

axis equal
set(gca,'visible','off')
colormap('jet')
quiver3([0 0 0],[0 0 0],[0 0 0],[1.3 0 0],[0 1.3 0],[0 0 1.3],'k-','LineWidth',2)
view(10,10)
xlim([-1 1.3])
ylim([-1 1.3])
zlim([zsurf 1.3])

% alpha 0.8

exportgraphics(gcf,fullfile(saveFolder,"tiltDRP_schematics.tif"),'Resolution',600)

%% exposure time profiles
figure("Units","centimeters","Position",[2 2 16 8])
tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
nexttile(1)
DRPdisp_intensity(bgDRP_00,exp_para_plot,cRange=[0 10],colormap="sky")
% colorbar("FontSize",14)
nexttile(2)
DRPdisp_intensity(bgDRP_40,exp_para_plot,cRange=[0 10],colormap="sky")
% colorbar("FontSize",14)
% nexttile(3)
% DRPdisp_intensity(bgDRP_00 ./ bgDRP_20 ,exp_para,cRange=[0 2])
% colorbar("FontSize",14)
exportgraphics(gcf,fullfile(saveFolder,"expTime_profile.tif"),'Resolution',300)

%% plot JET colormap colorbar
close all
plot_colorbar(cmap="jet",size=[750 50])
exportgraphics(gcf,fullfile(saveFolder,"jet_colorbar.tif"),'Resolution',300)
%% plot parula colormap colorbar
close all
plot_colorbar(cmap="parula",size=[750 50])
exportgraphics(gcf,fullfile(saveFolder,"parula_colorbar.tif"),'Resolution',300)

%% sample DRPs (measured, background, ground truth)
idx_xx = 70;
idx_yy = 300;
figure("Units","centimeters","Position",[2 2 18 6])
tiledlayout(1,3,'TileSpacing','tight','Padding','tight')
drpSample = splitDRP(exp_para,squeeze(igray_sample(idx_xx,idx_yy,:)),phitheta);
drpBack = splitDRP(exp_para,squeeze(igray_back(idx_xx,idx_yy,:)),phitheta);
drpNorm = uint8(double(drpSample)./double(drpBack)/2*256);
drpNorm(drpBack < 10) = 0;
nexttile(1)
DRPdisp(drpSample,exp_para)
nexttile(2)
DRPdisp(drpBack,exp_para)
nexttile(3)
DRPdisp(drpNorm,exp_para)
exportgraphics(gcf,fullfile(saveFolder,"measuredDRPs.tif"),'Resolution',300)
%% --------------------------------------------------------------------------
%% plot indexing performance of DI approach on full and partial DRPs
% plot maksed DRP and its reconstruction outcome

% drpLib = DRPLibGenerator(3*degree,exp_para);
thetaRange = 0:5:60;
outcome = cell(length(thetaRange),2);
phi = repmat(linspace(exp_para.ph_min, exp_para.ph_max, exp_para.ph_num),exp_para.th_num,1);
theta = repmat(linspace(exp_para.th_min, exp_para.th_max, exp_para.th_num)',1,exp_para.ph_num);
ii = 0;
for tiltAngle = thetaRange
    % mask = createTiltMask(tiltAngle,exp_para);
    ii = ii+1;
    mask = createTiltMask(tiltAngle,exp_para);
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

%% --------------------------------------------------------------------------
%% investigation of indexing performance on tilted samples
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
misOriPixel = cell2mat(misOriGrain_data);

%% area fraction calculation
tiltangle = 0:0.2:90;
n_tangle = length(tiltangle);
areaFrac = zeros(2,n_tangle);
HW = 8;

% generate datapoints on the unit sphere
resolution = 0.1;
S2G = equispacedS2Grid('resolution',resolution*degree,'maxTheta',90*degree);
phiGrid = (S2G.rho ./ degree)';
thetaGrid = (S2G.theta ./ degree)';
rGrid = [sind(thetaGrid).*cosd(phiGrid), sind(thetaGrid).*sind(phiGrid), cosd(thetaGrid)];

% generate [001] facet and reflectance peak
r_rot = ori_1deg * Miller({0,0,1},cs);
th_fz = r_rot.theta;
ph_fz = r_rot.rho;
th_fz_r = 2*th_fz;
x_fz = sin(th_fz_r).*cos(ph_fz);
y_fz = sin(th_fz_r).*sin(ph_fz);
z_fz = cos(th_fz_r);
pos_xyz = unique([x_fz, y_fz, z_fz],'row');
x_fz = pos_xyz(:,1);
y_fz = pos_xyz(:,2);
z_fz = pos_xyz(:,3);
notApproachable = pos_xyz(:,3) < sin((8+5-HW)*degree) | pos_xyz(:,3) > sin((65+HW)*degree);


% DRM setup mask
idx_DRMValid = thetaGrid >= 25 & thetaGrid <= 82;

for ii = 1:n_tangle
    % generate tilting surface normal direction
    alpha = tiltangle(ii);
    beta = 0;
    rN = [-sind(alpha)*cosd(beta), -sind(alpha)*sind(beta), cosd(alpha)];
    
    % calculate angle between rGrid and rN
    angleDev = acosd(rGrid * rN');
    angleDev2 = acosd(pos_xyz * rN');
    areaFrac(1,ii) = sum(angleDev(idx_DRMValid) > 90+HW);
    areaFrac(2,ii) = sum(angleDev2(~notApproachable) > 90+HW);
end
figure, plot(tiltangle, areaFrac(1,:)/length(phiGrid(idx_DRMValid)), ...
    tiltangle,areaFrac(2,:)/length(pos_xyz(~notApproachable,:)))

%% clear data00 data_temp misOriGrain_data
nn = 22892;
ori_ref = orientation.byEuler(drpLib.eulerDic.*degree,cs);
numPoints = size(outcome,1);
angleDiffTilt = zeros(nn,numPoints);
for ii = 1:numPoints
    ori_temp = orientation.byEuler(outcome{ii,1}.*degree,cs);
    angleDiffTilt(:,ii) = angle(ori_ref,ori_temp,cs)./degree;
end
clear ori_temp
figure("Units","centimeters","Position",[2 2 18 9])
plot(1:9,sum(angleDiffTilt > 5)/nn,'-','LineWidth',2,'Color','red')
hold on
set(gca,'LineWidth',2,'FontSize',14)
xticks(0:9)
xticklabels(0:5:45)
xlim([0 9.1])
ylim([0 0.6])
yticks(0:0.1:0.6)
yticklabels(0:10:60)
angleDiffTilt_111 = zeros(nn,numPoints);
% angleDiffTilt_110 = zeros(nn,numPoints);
for ii = 1:numPoints
    ori_temp = orientation.byEuler(outcome_111{ii,1}.*degree,cs);
    angleDiffTilt_111(:,ii) = angle(ori_ref,ori_temp,cs)./degree;
    % ori_temp = orientation.byEuler(outcome_110{ii,1}.*degree,cs);
    % angleDiffTilt_110(:,ii) = angle(ori_ref,ori_temp,cs)./degree;
end
plot(1:9,sum(angleDiffTilt_111 > 5)/nn,'-','LineWidth',2,'Color','blue')
% plot(1:9,sum(angleDiffTilt_110 > 5)/nn,':','Marker','square','LineWidth',2,'Color','black')

plot(0:9,1-sum(misOriPixel < 20)/size(misOriPixel,1)/0.9, ...
    ':','Marker','o','LineWidth',2,'Color','black')

% xxxx = (0:0.2:45)/5;
% plot(xxxx,areaFrac(2,1:length(xxxx))/length(pos_xyz(~notApproachable,:)),'-','LineWidth',2,"color",'red')

legend("Loss_{(100)} (SIM)","Loss_{(111)} (SIM)","Loss (EXP)",'Location','northwest')
xlabel("Tilt angle (deg)",'FontSize',14)
ylabel("Indexing Loss (%)",'FontSize',14)
exportgraphics(gcf,fullfile(saveFolder,"indexingSimPerformance_rr.tif"),'Resolution',600)
% print(gcf,fullfile(saveFolder,"indexingSimPerformance_R"),'-r600','-dtiffn')

%% experimental indexing error as a function of tilting angles
a = prctile(misOriPixel,25);
b = prctile(misOriPixel,50);
c = prctile(misOriPixel,75);
figure('unit','centimeters','position',[2 2 12 6])
% plot(x,a,'LineWidth',2,'Color','#fdae61')
hold on
plot(degrees,b,'x-','LineWidth',2,'Color','#3182bd')
% plot(x,c,'LineWidth',2,'Color','#fdae61')
fill([degrees, fliplr(degrees)], [a, fliplr(c)], [158,202,225]/255, 'FaceAlpha', 0.3, 'EdgeColor', [158,202,225]/255);
xlim([0 45.2])
ylim([0 62])
yticks(0:20:60)
xticks(0:5:45)
xlabel("Tilt angle (deg)")
ylabel("Indexing error (deg)")
box on
set(gca,'LineWidth',2,'FontSize',14)
% exportgraphics(gcf,fullfile(saveFolder,"indexingExpPerformance.tif"),'Resolution',600)
% print(gcf,fullfile(saveFolder,"indexingExpPerformance"),'-r600','-dtiffn')

%% make it a quasi-contour plot
% colorSeq = {'#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5',...
%     '#08519c','#08306b'};
colorSeq = [[247,251,255];[222,235,247];[198,219,239];[158,202,225];[107,174,214];...
    [66,146,198];[33,113,181];[8,81,156];[8,48,107]]./256;
colorSep = flipud(colorSeq);

lineColorSeq = [[254,229,217];[252,174,145];[251,106,74];[222,45,38];[165,15,21]]./256;
% lineColorSeq = flipud(lineColorSeq);
figure('unit','centimeters','position',[2 2 12 6])
hold on
layer_cur = zeros(1,10);
for ii = 1:9
    layer_nxt = prctile(misOriPixel,ii*10);
    % if mod(ii,2)==1
    %     plot(degrees,layer_nxt,'-','LineWidth',1,'Color','r')
    % end
    fill([degrees,fliplr(degrees)],[layer_cur,fliplr(layer_nxt)],colorSeq(ii,:),"EdgeColor",'none')
    layer_cur = layer_nxt;
end
for ii = 1:9
    layer_nxt = prctile(misOriPixel,ii*10);
    if mod(ii,2)==1
        plot(degrees,layer_nxt,'-','LineWidth',1,'Color',lineColorSeq((ii+1)/2,:))
    end
end
xlim([0 45.2])
ylim([0 62])
yticks(0:20:60)
xticks(0:5:45)
xlabel("Tilt angle (deg)")
ylabel("Indexing error (deg)")
box on
set(gca,'LineWidth',2,'FontSize',14)

% exportgraphics(gcf,fullfile(saveFolder,"errorHistoContour.tif"),Resolution=300)
%% --------------------------------------------------------------------------
%% tilt experiment errorMap
deg_temp = 0:10:40;
for ii = 1:5
    data_temp = load(fullfile("/Users/chenyangzhu/Desktop/datasetTemp", ...
        [num2str(deg_temp(ii),'%02d'),'deg'],"/indexed.mat"),"misOriMap");
    misOriMapTemp = data_temp.misOriMap;
    figure, imshow(misOriMapTemp,Border="tight")
    colormap(sky)
    clim([0 15])
    exportgraphics(gcf,fullfile(saveFolder,[num2str(deg_temp(ii),'%02d'),'deg_idxError.tif']),'Resolution',300)
    close all

end
clear deg_temp

%% grain categorizing and error heatmap
% degrees = 0:5:45;
% num_deg = length(degrees);
% misOriGrain_data = cell(1,num_deg);
% for ii = 1:num_deg
%     data_temp = load(fullfile("/Users/chenyangzhu/Desktop/datasetTemp", ...
%         [num2str(degrees(ii),'%02d'),'deg'],"/indexed.mat"),"misOriGrain");
%     data00 = data_temp.misOriGrain;
%     misOriGrain_data{ii} = data00(~isnan(data00));
%     workbar(ii/num_deg,sprintf("load %02d / %02d folders",[ii num_deg]));
% end
% clear data00 data_temp
% misOriGrains = cell2mat(misOriGrain_data);

misOriGrainAve = mean(misOriGrains,2);
misOriGrainStd = std(misOriGrains,1,2);
grainFeature = [misOriGrainAve, misOriGrainStd, max(misOriGrains,[],2),min(misOriGrains,[],2)];

grainTotOff = grainFeature(:,4) > 20;
grainWell = grainFeature(:,3) < 20;
grainJump = grainFeature(:,2) > 16;
grainStable = grainFeature(:,2) < 3;


figure("Units","centimeters","Position",[2 2 5*2 4*sum(grainJump)/10/2])
tiledlayout(1,2)
largeSTDGrains = misOriGrains(grainJump,:);
nSTD = size(largeSTDGrains,1);
idx_STD = find(grainJump);
nexttile(1)

s1 = pcolor(flipud(largeSTDGrains(1:floor(nSTD/2),:)));
s1.EdgeColor = 'none';
clim([0 62])
axis equal
xlim([1 10])
yticks(1:floor(nSTD/2))
yticklabels(flipud(idx_STD(1:floor(nSTD/2))))
xticks([])
colormap(sky)

nexttile(2)
s2 = pcolor(flipud(largeSTDGrains((floor(nSTD/2)+1):(nSTD-1),:)));
s2.EdgeColor = 'none';
clim([0 62])
axis equal
xlim([1 10])
yticks(1:floor(nSTD/2))
yticklabels(flipud(idx_STD((floor(nSTD/2)+1):(nSTD-1))))
xticks([])
colormap(sky)

exportgraphics(gcf,fullfile(saveFolder,"largeSTDHeatmap.tif"),'Resolution',300)

%% EBSD ground truth of Al sample
figure, imshow(plot_ipf_map(EUmap_ebsd),'Border','tight')
exportgraphics(gcf,fullfile(saveFolder,"EBSD_IPFZ.tif"),'Resolution',300)

%% Rotating back to local coordinates
euMap_ori = indexResult.euMap(:,pos_l:pos_r,:);
pos_l = 11;
pos_r = 850;
num_pxl = pos_r-pos_l+1;
rotAngle_list = linspace(45,-45,num_pxl);
euMap_temp = zeros(3,size(indexResult.euMap,1),num_pxl);
for ii = 1:num_pxl
    eu_list = squeeze(indexResult.euMap(:,pos_l+ii-1,:));
    ori_list = orientation.byEuler(eu_list*degree,cs);
    rot = rotation.byAxisAngle(vector3d.X,rotAngle_list(ii)*degree);
    ori_list_rot = rot * ori_list;
    euMap_temp(:,:,ii) = [ori_list_rot.phi1;ori_list_rot.Phi;ori_list_rot.phi2]./degree;
end
euMap_temp = permute(euMap_temp,[2,3,1]);
figure, imshowpair(plot_ipf_map(euMap_ori),plot_ipf_map(euMap_temp),'montage')

%% plot affected DRPs and crystal structure
exp_para_100 = exp_para;
exp_para_100.faceting = [1 0 0];
drpLib_100 = DRPLibGenerator(4*degree, exp_para_100);

tiltAngle = 30;
mask = createTiltMask(tiltAngle,exp_para);
% create masked DRP dictionary
drpDic_masked_100 = maskDRPLib(drpLib_100.drpDic,mask);

% perform indexing on masked simulated DRPs
nn = length(drpDic_masked_100);
euler = zeros(nn,3);
K = 3;
idx = zeros(nn,K);
batchNum = fix(nn/2);
[euler(1:batchNum,:),idx(1:batchNum,:)] = DirectDIEngine( ...
    drpDic_masked_100(1:batchNum), drpLib_100, K=K);
[euler(batchNum+1:nn,:),idx(batchNum+1:nn,:)] = DirectDIEngine( ...
    drpDic_masked_100(batchNum+1:nn), drpLib_100, K=K);
outcome_100{1,1} = euler;
outcome_100{1,2} = idx;

ori_ref_100 = orientation.byEuler(drpLib_100.eulerDic.*degree,cs);
ori_temp_100 = orientation.byEuler(outcome_100{1,1}.*degree,cs);
angleDiffTilt_100 = angle(ori_ref_100,ori_temp_100,cs)./degree;
% AE training
% drpDic_111_r = cell(1,length(drpLib_111.drpDic));
% for ii = 1:length(drpLib_111.drpDic)
%     drpDic_111_r{ii} = double(drpLib_111.drpDic{ii})/256;
% end
% hiddenSize1 = 100;
% AE_DRM = trainAutoencoder(drpDic_111_r,hiddenSize1, ...
%     'MaxEpochs',200, ...
%     'L2WeightRegularization',0.001, ...
%     'SparsityRegularization',4, ...
%     'SparsityProportion',0.10, ...
%     'ScaleData', false, ...
%     'UseGPU', false);
% clear hiddenSize1
% take out some of the grains with large indexing discrepancies
idx_list = find(angleDiffTilt_100 > 10);
idx = idx_list(randi(length(idx_list)));
% drp_recon = predict(AE_DRM,double(drpDic_masked{idx})/256);
figure('Units','centimeters','Position',[2 2 15 5])
tiledlayout(1,3,'TileSpacing','tight','Padding','tight')
nexttile(2)
DRPdisp(drpDic_masked_100{idx},exp_para)
% nexttile(2)
% DRPdisp(drp_recon,exp_para);
euler_temp = outcome_100{1};
nexttile(3)
DRPdisp(DRPsim(euler_temp(idx,1),euler_temp(idx,2),euler_temp(idx,3),exp_para_100),exp_para_100)
nexttile(1)
DRPdisp(drpLib_100.drpDic{idx},exp_para)
% nexttile(5)
% DRPdisp(predict(AE_DRM,double(drpLib_111.drpDic{idx})/256),exp_para)

exportgraphics(gcf,fullfile(saveFolder,"100_DRP_fail.tif"),Resolution=300)

% plot crystal structure
figure('Units','centimeters',Position=[2 2 5 5])
hold on
quiver3([0 0 0],[0 0 0],[0 0 0],[1.3 0 0],[0 1.3 0],[0 0 1.3],'k-','LineWidth',2)

% Cube (for {100} faceting)
r = 1;
V1=[-1;  1; 1; -1; -1;  1; 1; -1;];
V2=[-1; -1; 1;  1; -1; -1; 1;  1;];
V3=[-1; -1;-1; -1;  1;  1; 1;  1;];
F= [1 2 3 4; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 5 6 7 8;];
[THETA,PHI,R]=cart2sph(V1,V2,V3);
R=r.*ones(size(V1(:,1)));
[V1,V2,V3]=sph2cart(THETA,PHI,R);
V=[V1 V2 V3];
% apply the rotation
euTemp = drpLib_100.eulerDic(idx,:);
eu1 = euTemp(1);
eu2 = euTemp(2);
eu3 = euTemp(3);
roll = eu1*degree; pitch = eu2*degree; yaw = eu3*degree;
dcm = angle2dcm(roll, pitch, yaw, 'ZXZ');
V = V*dcm;
patch('Faces',F,'Vertices',V,'FaceColor',[227,26,28]/255,'FaceAlpha',0.818,'EdgeColor',[128,0,38]/255,...
    'LineWidth',2); axis equal; grid on; hold on; view(3); axis off;
xticks([])
yticks([])
zticks([])
xlim([-1,1])
ylim([-1,1])
zlim([-1,1])
view([37.5,30])
exportgraphics(gcf,fullfile(saveFolder,"100_original.tif"),Resolution=300)
%% --------------------------------------------------------------------------
%% functions used in this script
function [bgDRP] = createBG(el_angle,exp_para)
    th_max = exp_para.th_max;
    th_min = exp_para.th_min;
    th_num = exp_para.th_num;
    ph_num = exp_para.ph_num;
    th_step = (th_max - th_min) / (th_num - 1);
    ph_step = 360 / ph_num;

    th_range = th_min : th_step : th_max;
    th_DRP=repmat(transpose(th_range),1,ph_num);
    ph_DRP=repmat(0:ph_step:360-ph_step,length(th_range),1);

    vec_DRP=zeros(3,length(th_range),ph_num);
    bgDRP = zeros(th_num,ph_num);
    vecNorm = [sind(el_angle),0,cosd(el_angle)];
    
    for ii = 1:th_num
        for jj = 1:ph_num
            tmp_vec = thph2vec(th_DRP(ii,jj), ph_DRP(ii,jj));
            vec_DRP(:,ii,jj) = normr(tmp_vec);
            tmp_ang = acosd(vecNorm * normr(tmp_vec)');
            bgDRP(ii,jj) = sind(45+el_angle) / sind(90-tmp_ang);
            % bgDRP(ii,jj) = vecNorm * normr(tmp_vec)';
        end
    end  
end

function DRPdisp_intensity(drp,exp_para,options)
% this function is to display a drp in polar coordinates
% drp is in size of th_num x ph_num
% Create date: Aug 27, 2021
% Edit date: Sep 14, 2021
% By: Chenyang ZHU @ NTU
% -------------------------------------------------------------------------
arguments
    drp {mustBeNumeric}
    exp_para struct
    options.EdgeColor (1,1) string = 'none'
    options.colormap (1,1) string = 'jet'
    options.project (1,1) string = "stereo"
    options.cRange (1,2) double = [0 255]
    options.scaleBar (1,1) logical = 0
%     options.title (1,1) string = ""
end

th_max = exp_para.th_max;
th_min = exp_para.th_min;
th_num = exp_para.th_num;
ph_num = exp_para.ph_num;

th_step = (th_max - th_min) / (th_num - 1);
ph_step = 360 / (ph_num - 1);

[x,y] = meshgrid(0:ph_step:360,th_min:th_step:th_max);
if options.project == "stereo"
    xx = cosd(y).*cosd(x)./(1+sind(y));
    yy = cosd(y).*sind(x)./(1+sind(y));
elseif options.project == "direct"
    xx = cosd(y).*cosd(x);
    yy = cosd(y).*sind(x);
end
h = pcolor(xx,yy,drp);
axis tight
axis equal
set(h, 'EdgeColor', options.EdgeColor);
ax = gca;
set(ax,'visible','off')
colormap(ax,options.colormap)
% title(h,options.title,'FontSize',16)
caxis(options.cRange)
if options.scaleBar
    colorbar('FontSize',16)
end
end

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

